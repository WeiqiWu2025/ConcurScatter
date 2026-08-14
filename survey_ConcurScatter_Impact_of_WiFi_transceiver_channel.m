clc;
clear;
close all;
addpath(genpath(pwd));
rng(1993);

%%%%%*** Waveform Configuration ***%%%%%
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20';
cfgHT.NumTransmitAntennas = 1;
cfgHT.NumSpaceTimeStreams = 1;
cfgHT.PSDULength = 500;
cfgHT.MCS = 0;
cfgHT.ChannelCoding = 'BCC';

fs = wlanSampleRate(cfgHT); %#ok<NASGU>
ofdmInfo = wlanHTOFDMInfo('HT-Data',cfgHT); %#ok<NASGU>
ind = wlanFieldIndices(cfgHT);

tgnChannel_TX_RX = wlanTGnChannel;
tgnChannel_TX_RX.DelayProfile = 'Model-B';
tgnChannel_TX_RX.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgnChannel_TX_RX.NumReceiveAntennas = 1;
tgnChannel_TX_RX.TransmitReceiveDistance = 1;
tgnChannel_TX_RX.LargeScaleFadingEffect = 'None';

tgnChannel_Tag_RX = wlanTGnChannel;
tgnChannel_Tag_RX.DelayProfile = 'Model-B';
tgnChannel_Tag_RX.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgnChannel_Tag_RX.NumReceiveAntennas = 1;
tgnChannel_Tag_RX.TransmitReceiveDistance = 1;
tgnChannel_Tag_RX.LargeScaleFadingEffect = 'None';

%%%%%*** Simulation Parameters ***%%%%%
snr = 10:10:30;
global numTags;
numTags = 3;

tag_modulation = 2;
len_fft = 64;
len_cp = 16;
num_circshift = linspace(0,len_fft-1,numTags).';
num_circshift1 = floor(num_circshift);
num_circshift2 = num_circshift-num_circshift1;

requestedSeqLenForEstChannel = 20;
maxNumPackets = 30000;

S = numel(snr);
numBitErrs = zeros(S,numTags);
berEst = zeros(S,numTags);

WiFi_transceiver_link = 'high quality';
packetLoss = zeros(S,1);

for i = 1:S
    disp(['SNR: ',num2str(snr(i)),' dB...']);
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = 1;
    RandStream.setGlobalStream(stream);

    n = 1;
    while n<=maxNumPackets
        disp(['SNR: ',num2str(snr(i)),' dB -> ',...
            'n: ',num2str(n),'-th packet']);

        txPSDU = randi([0 1],cfgHT.PSDULength*8,1);
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx;zeros(15,cfgHT.NumTransmitAntennas)];

        exSig = [];
        for tagIdx = 1:numTags
            coefficient = 0.1*(-1+2*rand)+1i*0.1*(-1+2*rand);
            exSig = [exSig,tx.*coefficient]; %#ok<AGROW>
        end

        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        if mod(temp,2)==1
            numTagData = (numel(tx)-720-15-80-80-80)/80;
        else
            numTagData = (numel(tx)-720-15-80-80)/80;
        end

        [tagData,payloadTagData,seqLenForEstChannel,...
            ~,numPayload] = ...
            survey_ConcurScatter_funcPrepareCEFrame(...
            numTagData,numTags,requestedSeqLenForEstChannel,...
            tag_modulation);

        timeDomainTagData = repelem(tagData,1,len_cp+len_fft);
        fShift = (num_circshift+len_fft)./len_fft;
        sampleIndex = 1:len_cp+len_fft;
        shiftSignal = repmat(...
            exp(1i*2*pi*fShift*(sampleIndex-1)),1,numTagData);
        tagModulationSignal = shiftSignal.*timeDomainTagData;
        tagModulationSignal = [zeros(numTags,len_cp),...
            tagModulationSignal(:,1:end-len_cp)].';

        bxSig = cell(1,numTags);
        for tagIdx = 1:numTags
            bxSig{tagIdx} = exSig(:,tagIdx);
            bxSig{tagIdx}(801:800+length(tagModulationSignal)) = ...
                tagModulationSignal(:,tagIdx).*...
                bxSig{tagIdx}(801:800+length(tagModulationSignal));
        end

        for tagIdx = 1:numTags
            reset(tgnChannel_Tag_RX);
            bxSig{tagIdx} = tgnChannel_Tag_RX(bxSig{tagIdx});
        end

        reset(tgnChannel_TX_RX);
        rxSigFromWiFiTX = tgnChannel_TX_RX(tx);
        WiFiLinkSNR = survey_ConcurScatter_funcChannelQuality2SNR(...
            WiFi_transceiver_link);

        rx = complex(zeros(length(bxSig{1}),1));
        for tagIdx = 1:numTags
            rx = rx+bxSig{tagIdx};
        end

        [rxFromTags,~,~] = func_awgn(rx,snr(i),'measured');
        [rxFromWiFiTX,~,~] = func_awgn(...
            rxSigFromWiFiTX,WiFiLinkSNR,'measured');
        [decodedPSDU,flag] = ...
            survey_ConcurScatter_funcWiFiRX(rxFromWiFiTX,cfgHT);

        if flag==0
            packetLoss(i) = packetLoss(i)+1;
            % Only payload bits are part of the reported BER. Reference
            % symbols carry no payload and are not counted as bit errors.
            numBitErrs(i,:) = numBitErrs(i,:)+numPayload;
            n = n+1;
            continue;
        end

        ofdmDemod = survey_ConcurScatter_funcReceiver(...
            rxFromTags(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        receivedTagSymbols = ofdmDemod(:,2:1+numTagData);

        [~,originalSymbols] = ...
            survey_ConcurScatter_funcOFDMSymDerived(decodedPSDU,cfgHT);
        excitationSymbols = originalSymbols(:,2:1+numTagData);

        [decodedSymbols,H_est] = ... %#ok<NASGU>
            survey_ConcurScatter_funcChannelAwareDecode(...
            receivedTagSymbols,excitationSymbols,numTags,...
            seqLenForEstChannel,num_circshift1,num_circshift2);
        decodedBits = pskdemod(decodedSymbols,tag_modulation);
        originalBits = pskdemod(payloadTagData,tag_modulation).';

        for tagIdx = 1:numTags
            numBitErrs(i,tagIdx) = numBitErrs(i,tagIdx)+...
                biterr(originalBits(:,tagIdx),decodedBits(:,tagIdx));
        end
        n = n+1;
    end

    for tagIdx = 1:numTags
        berEst(i,tagIdx) = numBitErrs(i,tagIdx)/...
            (numPayload*maxNumPackets);
    end
end

aaa = 1; %#ok<NASGU>



