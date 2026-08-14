clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993); % For repeatable results

%%%%%*** Waveform Configuration ***%%%%%
% Create a format configuration object for a 1-by-1 HT transmission.
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1;
cfgHT.NumSpaceTimeStreams = 1;
cfgHT.PSDULength = 1000; % PSDU length in bytes
cfgHT.MCS = 0; % One spatial stream, BPSK, rate 1/2
cfgHT.ChannelCoding = 'BCC';

fs = wlanSampleRate(cfgHT); %#ok<NASGU> % Baseband sampling rate
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); %#ok<NASGU>
ind = wlanFieldIndices(cfgHT);

% Create and configure the tag-to-receiver channel.
tgnChannel = wlanTGnChannel;
tgnChannel.DelayProfile = 'Model-B';
tgnChannel.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgnChannel.NumReceiveAntennas = 1;
tgnChannel.TransmitReceiveDistance = 1;
tgnChannel.LargeScaleFadingEffect = 'None';

%%%%%*** Simulation Parameters ***%%%%%
snr = 30;
global numTags;
numTags = 2;

tag_modulation = 8;
len_fft = 64;
len_cp = 16;
num_circshift = linspace(0,len_fft-1,numTags).';
num_circshift1 = floor(num_circshift);
num_circshift2 = num_circshift-num_circshift1;

% Mecha-style tag-specific reference fields. During tag i's reference
% field, only tag i backscatters a known +1 symbol. The tag still applies
% its ConcurScatter cyclic shift, so the LS regressor is the reconstructed
% (shifted) pattern rather than the unshifted excitation symbol.
seqLenForEstChannel = 20;

maxNumPackets = 1000;

S = numel(snr);
numBitErrs = zeros(S,numTags);
berEst = zeros(S,numTags);

for i = 1:S
    disp(['SNR: ',num2str(snr(i)),' dB...']);

    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = 1;
    RandStream.setGlobalStream(stream);

    n = 1;
    while n<=maxNumPackets
        disp(['snr: ',num2str(snr(i)),' dB -> ','n: ',num2str(n),'-th packet']);

        %%%%%*** TX side ***%%%%%
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1);
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx;zeros(15,cfgHT.NumTransmitAntennas)];

        exSig = [];
        %%%%%*** TX-to-tag backscatter channel ***%%%%%
        for chan_tx_tag_idx1 = 1:numTags
            bxCoeffForTxTag_real = -1+2*rand(1,1);
            bxCoeffForTxTag_imag = -1+2*rand(1,1);
            bxCoeffForTxTag = 0.1*bxCoeffForTxTag_real+...
                1i*0.1*bxCoeffForTxTag_imag;
            tmp_exSig = tx.*bxCoeffForTxTag;
            exSig = [exSig,tmp_exSig]; %#ok<AGROW>
        end

        %%%%%*** Tags side ***%%%%%
        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        if mod(temp,2)==1
            numSymForPsdu = (numel(tx)-720-15-80-80-80)/80;
        else
            numSymForPsdu = (numel(tx)-720-15-80-80)/80;
        end
        numTagData = numSymForPsdu;

        numReferenceSymbols = seqLenForEstChannel*numTags;
        numPayload = numTagData-numReferenceSymbols;
        if numPayload<=0
            error(['The HT-Data field is too short for the requested ',...
                'tag-specific reference fields.']);
        end

        % tagData has dimensions numTags-by-numTagData. The first
        % seqLenForEstChannel*numTags symbols are orthogonal reference
        % fields; the remaining symbols are concurrent payload data.
        tagData = complex(zeros(numTags,numTagData));
        for tag_idx = 1:numTags
            refCols = (tag_idx-1)*seqLenForEstChannel+...
                (1:seqLenForEstChannel);
            tagData(tag_idx,refCols) = 1;
        end
        payloadTagData = reshape(...
            survey_ConcurScatter_funcRandd(...
            numTags*numPayload,tag_modulation),numPayload,[]).';
        tagData(:,numReferenceSymbols+1:end) = payloadTagData;

        time_domain_tagData = repelem(tagData,1,len_cp+len_fft);
        f_shift = (num_circshift+len_fft)./len_fft;
        tmp_time_domain_samples = 1:len_cp+len_fft;
        sig_shift = repmat(...
            exp(1i*2*pi*f_shift*(tmp_time_domain_samples-1)),...
            1,numTagData);
        sig_mod = sig_shift.*time_domain_tagData;
        sig_mod = [zeros(numTags,len_cp),sig_mod(:,1:end-len_cp)];
        sig_mod = sig_mod.';

        bxSig = cell(1,numTags);
        for tag_idx1 = 1:numTags
            bxSig{tag_idx1} = exSig(:,tag_idx1);
            bxSig{tag_idx1}(801:800+length(sig_mod)) = ...
                sig_mod(:,tag_idx1).*...
                bxSig{tag_idx1}(801:800+length(sig_mod));
        end

        %%%%%*** Tag-to-receiver backscatter channel ***%%%%%
        for chan_tag_rx_idx1 = 1:numTags
            reset(tgnChannel);
            bxSig{chan_tag_rx_idx1} = ...
                tgnChannel(bxSig{chan_tag_rx_idx1});
        end

        %%%%%*** RX side ***%%%%%
        rx = complex(zeros(length(bxSig{1}),1));
        for rx_idx1 = 1:numTags
            rx = rx+bxSig{rx_idx1};
        end

        [rxFromTags,~,~] = func_awgn(rx,snr(i),'measured');
        ofdmDemod = survey_ConcurScatter_funcReceiver(...
            rxFromTags(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        ofdmDemod1 = ofdmDemod(:,2:1+numTagData);

        [~,original_ofdm_symbols] = ...
            survey_ConcurScatter_funcOFDMSymDerived(txPSDU,cfgHT);
        original_ofdm_symbols1 = ...
            original_ofdm_symbols(:,2:1+numTagData);

        %%%%%*** Per-tag, per-subcarrier LS channel estimation ***%%%%%
        % H_est(row,column) = H_est(subcarrier,tag). It is indexed by the
        % received FFT-bin order. It must not be cyclically shifted again.
        H_est = complex(zeros(len_fft,numTags));

        for tag_idx = 1:numTags
            refCols = (tag_idx-1)*seqLenForEstChannel+...
                (1:seqLenForEstChannel);

            receivedReference = ofdmDemod1(:,refCols);
            excitationReference = original_ofdm_symbols1(:,refCols);

            pos = num_circshift1(tag_idx);
            rem = num_circshift2(tag_idx);
            referencePattern = reconstructConcurScatterPattern(...
                excitationReference,pos,rem);

            numerator = sum(...
                conj(referencePattern).*receivedReference,2);
            denominator = sum(abs(referencePattern).^2,2);
            valid = denominator>1e-12;

            H_est(valid,tag_idx) = ...
                numerator(valid)./denominator(valid);
        end

        %%%%%*** Channel-aware concurrent payload decoding ***%%%%%
        payloadCols = numReferenceSymbols+1:numTagData;
        receivedPayload = ofdmDemod1(:,payloadCols);
        excitationPayload = original_ofdm_symbols1(:,payloadCols);

        dem_rec = complex(zeros(numTags,numPayload));
        org_circ = complex(zeros(numTags,len_fft,numPayload));

        for rx_idx2 = 1:numTags
            pos = num_circshift1(rx_idx2);
            rem = num_circshift2(rx_idx2);

            idealPattern = reconstructConcurScatterPattern(...
                excitationPayload,pos,rem);

            % The cyclic shift is already contained in idealPattern.
            % H_est remains in received FFT-bin order.
            channelAwarePattern = ...
                H_est(:,rx_idx2).*idealPattern;

            org_circ(rx_idx2,:,:) = channelAwarePattern;
            dem_rec(rx_idx2,:) = sum(...
                receivedPayload.*conj(channelAwarePattern),1);
        end

        % Construct C_H(k,a)=<Q_a,Q_k>, where
        % Q_a=H_est(:,a).*P_a is the channel-aware pattern of tag a.
        conj_org = complex(zeros(numTags,numTags,numPayload));
        for rx_idx3 = 1:numTags
            for rx_idx4 = 1:numTags
                conj_org(rx_idx3,rx_idx4,:) = sum(...
                    org_circ(rx_idx4,:,:).*...
                    conj(org_circ(rx_idx3,:,:)),2);
            end
        end

        % Solve A=C_H*d independently for every payload OFDM symbol.
        dec_tag = complex(zeros(numPayload,numTags));
        for payload_idx = 1:numPayload
            correlationVector = dem_rec(:,payload_idx);
            correlationMatrix = ...
                reshape(conj_org(:,:,payload_idx),numTags,numTags);
            decodedSymbols = pinv(correlationMatrix)*correlationVector;
            dec_tag(payload_idx,:) = decodedSymbols.';
        end

        dec_tag = pskdemod(dec_tag,tag_modulation);

        org_tag = pskdemod(payloadTagData,tag_modulation).';

        % Calculate payload BER; reference symbols are not counted.
        for tt1 = 1:numTags
            numBitErrs(i,tt1) = numBitErrs(i,tt1)+...
                biterr(org_tag(:,tt1),dec_tag(:,tt1));
        end
        n = n+1;
    end

    for tt2 = 1:numTags
        berEst(i,tt2) = numBitErrs(i,tt2)/...
            (numPayload*maxNumPackets);
    end
end

aaa = 1; %#ok<NASGU>


function pattern = reconstructConcurScatterPattern(excitationSymbols,pos,rem)
% Reconstruct the integer/fractional cyclic-shift pattern used by the
% original ConcurScatter implementation. Rows are received FFT bins and
% columns are OFDM symbols.

pattern = ...
    circshift(excitationSymbols,pos)*...
        sinc(-rem)*exp(-1i*pi*-rem)+...
    circshift(excitationSymbols,pos+1)*...
        sinc(1-rem)*exp(-1i*pi*(1-rem))+...
    circshift(excitationSymbols,pos+2)*...
        sinc(2-rem)*exp(-1i*pi*(2-rem))+...
    circshift(excitationSymbols,pos-1)*...
        sinc(-1-rem)*exp(-1i*pi*(-1-rem))+...
    circshift(excitationSymbols,pos+3)*...
        sinc(3-rem)*exp(-1i*pi*(3-rem))+...
    circshift(excitationSymbols,pos-2)*...
        sinc(-2-rem)*exp(-1i*pi*(-2-rem))+...
    circshift(excitationSymbols,pos+4)*...
        sinc(4-rem)*exp(-1i*pi*(4-rem))+...
    circshift(excitationSymbols,pos-3)*...
        sinc(-3-rem)*exp(-1i*pi*(-3-rem));
end
