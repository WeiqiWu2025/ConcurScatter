clc; 
clear; 
close all;
addpath(genpath(pwd)); 
rng(1993); % For repeatable results


%%%%%*** Waveform Configuration ***%%%%%
% Create a format configuration object for a 1-by-1 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 1; % 1 transmit antennas
cfgHT.NumSpaceTimeStreams = 1; % 1 space-time streams
cfgHT.PSDULength = 2000; % PSDU length in bytes % 64 %   PSUDLength * 8/26
cfgHT.MCS = 0; % 1 spatial streams, BPSK rate-1/2
cfgHT.ChannelCoding = 'BCC'; % BCC channel coding

fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); % Get the OFDM info
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet

%%%%%*** Simulation Parameters ***%%%%%
snr = 10:10:30; %10:2:30;
global numTags;
numTags = 3;

tag_modulation = 2; % 标签采用的PSK调制阶数：2=BPSK
len_fft = 64;
len_cp = 16;
num_circshift = (linspace(0,len_fft-1,numTags).');
num_circshift1 = floor(num_circshift); % % 提取CS的整数部分（CS_i）：对应循环移位向下取整
num_circshift2 = num_circshift - num_circshift1; % % 提取CS的分数部分（CS_f）：总CS-整数CS_i


maxNumPackets = 2000; % The maximum number of packets at an SNR point

S = numel(snr); % 返回数组snr中元素的个数
numBitErrs = zeros(S,numTags);
berEst = zeros(S,numTags);


for i = 1:S
    disp(['SNR: ',num2str(snr(i)),' dB...']);
    % Set random substream index per iteration to ensure that each
    % iteration uses a repeatable set of random numbers
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = i;
    RandStream.setGlobalStream(stream);
    
    % Loop to simulate multiple packets
    n = 1; % Index of packet transmitted
    while  n<=maxNumPackets
        disp(['snr: ',num2str(snr(i)),' dB -> ','n: ',num2str(n),'-th packet']);
        %%%%%*** TX side ***%%%%%
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)]; % Add trailing zeros to allow for channel filter delay
        
        exSig = [];
        H_TX_Tags = [];
        %%%%%*** TX-Tags backscatter channel & AWGN
        for chan_tx_tag_idx1 = 1:numTags
            bxCoeffForTxTag_real = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTxTag_real = bxCoeffForTxTag_real*0.1;
            bxCoeffForTxTag_imag = bxCoeffForTxTag_imag*0.1;
            bxCoeffForTxTag = bxCoeffForTxTag_real+1i*bxCoeffForTxTag_imag;
            tmp_exSig = tx.*bxCoeffForTxTag;
            exSig = [exSig,tmp_exSig];
            H_TX_Tags = [H_TX_Tags,bxCoeffForTxTag];
        end
        
        %%%%%*** Tags side ***%%%%%
        % Backscatter at the tag
        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        numSymForPsdu = 0;
        numSymForTailPad = 0;
        if mod(temp,2) == 1
            numSymForPsdu = (numel(tx)-720-15-80-80-80)/80;
            numSymForTailPad = 2;
        else
            numSymForPsdu = (numel(tx)-720-15-80-80)/80;
            numSymForTailPad = 1;
        end
        numTagData = numSymForPsdu; % modulate one tag data per one symbol
        
        % Initial tags data
        tagData = reshape(survey_ConcurScatter_funcRandd(numTags*numTagData,tag_modulation),numTagData,[]).';
        time_domain_tagData = repelem(tagData,1,len_cp+len_fft);
        f_shift = (num_circshift+len_fft)./len_fft;
        tmp_time_domain_samples = 1:1:len_cp+len_fft;
        sig_shift = repmat(exp(1i*2*pi*f_shift*(tmp_time_domain_samples-1)),1,numTagData);
        sig_mod = sig_shift.*time_domain_tagData;
        sig_mod = [zeros(numTags,len_cp),sig_mod(:,1:end-len_cp)];
        sig_mod = sig_mod.';
        
        for tag_idx1 = 1:numTags
            bxSig{tag_idx1} = exSig(:,tag_idx1);
            bxSig{tag_idx1}(801:800+length(sig_mod)) = sig_mod(:,tag_idx1).*bxSig{tag_idx1}(801:800+length(sig_mod));
        end
        
        %%%%%***** Backscatter channel & AWGN ***%%%%%
        H_Tags_RX = [];
        for chan_tag_rx_idx1 = 1:numTags
            bxCoeffForTagRx_real = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_imag = -1+(1+1)*rand(1,1);
            bxCoeffForTagRx_real = bxCoeffForTagRx_real*0.01;
            bxCoeffForTagRx_imag = bxCoeffForTagRx_imag*0.01;
            bxCoeffForTagRx = bxCoeffForTagRx_real+1i*bxCoeffForTagRx_imag;
            bxSig{chan_tag_rx_idx1} = bxSig{chan_tag_rx_idx1}.*bxCoeffForTagRx; % backscatter channel
            H_Tags_RX = [H_Tags_RX,bxCoeffForTagRx];
        end
        actualH = H_TX_Tags.*H_Tags_RX;
            
        %%%%%*** RX side ***%%%%%
        rx = complex(zeros(length(bxSig{1}),1));
        for rx_idx1 = 1:numTags
            rx = rx + bxSig{rx_idx1};
        end
%         rxFromTags = rx;
        [rxFromTags,~,~] = func_awgn(rx,snr(i),'measured'); % Received signal from Tags-RX channel
        ofdmDemod = survey_ConcurScatter_funcReceiver(rxFromTags(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        ofdmDemod1 = ofdmDemod(:,2:1+numTagData);
        
        [~,original_ofdm_symbols] = survey_ConcurScatter_funcOFDMSymDerived(txPSDU,cfgHT);
        original_ofdm_symbols1 = original_ofdm_symbols(:,2:1+numTagData);
        
        dem_rec = [];
        org_circ = [];
        for rx_idx2 = 1:numTags
            pos = num_circshift1(rx_idx2);
            rem = num_circshift2(rx_idx2);
            temp = circshift(original_ofdm_symbols1,pos) * sinc(-rem)*exp(-1i*pi*-rem) + ... 
                circshift(original_ofdm_symbols1,pos+1) * sinc(1-rem)*exp(-1i*pi*(1-rem)) + ... 
                circshift(original_ofdm_symbols1,pos + 2 ) * sinc(2-rem)*exp(-1i*pi*(2-rem))+... 
                circshift(original_ofdm_symbols1,pos -1 ) * sinc(-1-rem)*exp(-1i*pi*(-1-rem))+... 
                circshift(original_ofdm_symbols1,pos + 3 ) * sinc(3-rem)*exp(-1i*pi*(3-rem))+... 
                circshift(original_ofdm_symbols1,pos -2 ) * sinc(-2-rem)*exp(-1i*pi*(-2-rem))+... 
                circshift(original_ofdm_symbols1,pos + 4 ) * sinc(4-rem)*exp(-1i*pi*(4-rem))+... 
                circshift(original_ofdm_symbols1,pos - 3 ) * sinc(-3-rem)*exp(-1i*pi*(-3-rem)); 
            org_circ(rx_idx2,:,:) = temp;
            dem_rec(rx_idx2,:) = sum(ofdmDemod1.*conj(temp));
        end
        
        conj_org = [];
        for rx_idx3 = 1:numTags
            for rx_idx4 = 1:numTags
                conj_org(rx_idx3,rx_idx4,:) = sum(org_circ(rx_idx4,:,:).*conj(org_circ(rx_idx3,:,:)),2);
            end
        end
        
        matrix_A = [];
        temp = [];
        temp1 = [];
        for rx_idx5 = 1:numTags
            temp(:,:) = conj_org(rx_idx5,:,:); 
            temp1(:,:) = [dem_rec(rx_idx5,:);temp]; 
            matrix_A = [matrix_A;temp1];
        end
        matrix_A = reshape(matrix_A,numTags+1,[]).';
        
        dec_tag = [];
        for rx_idx6 = 0:numTagData-1
            dec_tag(rx_idx6+1,:) =  pinv(matrix_A(rx_idx6*numTags+[1:numTags],2:numTags+1))*matrix_A(rx_idx6*numTags+[1:numTags],1);
        end
        
        dec_tag = pskdemod(dec_tag,tag_modulation);
        
        org_tag = pskdemod(tagData,tag_modulation);
        org_tag = org_tag.';
        
        % calculate the number of bits
        for tt1 = 1:numTags
            numBitErrs(i,tt1) = numBitErrs(i,tt1) + biterr(org_tag(:,tt1),dec_tag(:,tt1));
        end
        n = n+1;
        
    end
    % calculate bit error rate
    for tt2 = 1:numTags
        berEst(i,tt2) = numBitErrs(i,tt2)/(numTagData*maxNumPackets);
    end
    
end

aaa = 1;


