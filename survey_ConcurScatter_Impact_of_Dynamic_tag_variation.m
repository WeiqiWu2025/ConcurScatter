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

tgnChannel = wlanTGnChannel;
tgnChannel.DelayProfile = 'Model-B';
tgnChannel.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgnChannel.NumReceiveAntennas = 1;
tgnChannel.TransmitReceiveDistance = 1;
tgnChannel.LargeScaleFadingEffect = 'None';

%%%%%*** Simulation Parameters ***%%%%%
snr = 30;

% The system is initially configured for numTags tags. numAddTags tags
% subsequently join without changing the reference/payload timing of the
% initial tags. Change numAddTags to 0, 1, or 2 for the three cases used in
% the dynamic-tag experiment.
global numTags;
numTags = 2;
global numAddTags;
numAddTags = 2;
numTotalTags = numTags+numAddTags;

tag_modulation = 2;
len_fft = 64;
len_cp = 16;
num_circshift = linspace(0,len_fft-1,numTotalTags).';
num_circshift1 = floor(num_circshift);
num_circshift2 = num_circshift-num_circshift1;

seqLenForEstChannel = 20;
maxNumPackets = 1000;

S = numel(snr);
numBitErrs = zeros(S,numTotalTags);
berEst = zeros(S,numTotalTags);

for i = 1:S
    disp(['SNR: ',num2str(snr),' dB...']);
    stream = RandStream('combRecursive','Seed',0);
    stream.Substream = 1;
    RandStream.setGlobalStream(stream);

    n = 1;
    while n<=maxNumPackets
        disp(['snr: ',num2str(snr),' dB -> ',...
            'n: ',num2str(n),'-th packet']);

        txPSDU = randi([0 1],cfgHT.PSDULength*8,1);
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        tx = [tx;zeros(15,cfgHT.NumTransmitAntennas)];

        exSig = complex(zeros(numel(tx),numTotalTags));
        for tagIdx = 1:numTotalTags
            coefficient = 0.1*(-1+2*rand)+1i*0.1*(-1+2*rand);
            exSig(:,tagIdx) = tx.*coefficient;
        end

        temp = ceil((cfgHT.PSDULength*8+16+6)/26);
        if mod(temp,2)==1
            numTagData = (numel(tx)-720-15-80-80-80)/80;
        else
            numTagData = (numel(tx)-720-15-80-80)/80;
        end

        numInitialReferenceSymbols = seqLenForEstChannel*numTags;
        numExpandedReferenceSymbols = ...
            seqLenForEstChannel*numTotalTags;
        numInitialPayload = numTagData-numInitialReferenceSymbols;
        numCommonPayload = numTagData-numExpandedReferenceSymbols;
        numInitialPayloadSymbolsLost = ...
            numExpandedReferenceSymbols-numInitialReferenceSymbols;

        if numCommonPayload<1
            error('ConcurScatter:InsufficientDynamicTagFrame',...
                ['The HT-Data field must contain the expanded reference ',...
                'region and at least one common payload symbol.']);
        end

        % Initial tags retain the frame constructed for the original tag
        % set: their reference field ends after L*numTags symbols and their
        % payload begins immediately afterwards.
        tagData = complex(zeros(numTotalTags,numTagData));
        initialReference = complex(zeros(...
            numTags,numInitialReferenceSymbols));
        for tagIdx = 1:numTags
            refCols = (tagIdx-1)*seqLenForEstChannel+...
                (1:seqLenForEstChannel);
            initialReference(tagIdx,refCols) = 1;
        end
        initialPayloadTagData = reshape(...
            survey_ConcurScatter_funcRandd(...
            numTags*numInitialPayload,tag_modulation),...
            numInitialPayload,[]).';
        tagData(1:numTags,:) = ...
            [initialReference,initialPayloadTagData];

        % Added tags use their positions in the expanded reference layout.
        % During these new reference slots, the initial tags are already
        % transmitting payload. This interference is intentionally retained
        % to reproduce the fixed-reference-field experiment in Mecha.
        if numAddTags>0
            expandedReference = complex(zeros(...
                numTotalTags,numExpandedReferenceSymbols));
            for tagIdx = numTags+1:numTotalTags
                refCols = (tagIdx-1)*seqLenForEstChannel+...
                    (1:seqLenForEstChannel);
                expandedReference(tagIdx,refCols) = 1;
            end
            addedPayloadTagData = reshape(...
                survey_ConcurScatter_funcRandd(...
                numAddTags*numCommonPayload,tag_modulation),...
                numCommonPayload,[]).';
            tagData(numTags+1:numTotalTags,:) = ...
                [expandedReference(numTags+1:numTotalTags,:),...
                addedPayloadTagData];
        else
            addedPayloadTagData = complex(zeros(0,numCommonPayload));
        end

        timeDomainTagData = repelem(tagData,1,len_cp+len_fft);
        fShift = (num_circshift+len_fft)./len_fft;
        sampleIndex = 1:len_cp+len_fft;
        shiftSignal = repmat(...
            exp(1i*2*pi*fShift*(sampleIndex-1)),1,numTagData);
        tagModulationSignal = shiftSignal.*timeDomainTagData;
        tagModulationSignal = [zeros(numTotalTags,len_cp),...
            tagModulationSignal(:,1:end-len_cp)].';

        bxSig = cell(1,numTotalTags);
        for tagIdx = 1:numTotalTags
            bxSig{tagIdx} = exSig(:,tagIdx);
            bxSig{tagIdx}(801:800+length(tagModulationSignal)) = ...
                tagModulationSignal(:,tagIdx).*...
                bxSig{tagIdx}(801:800+length(tagModulationSignal));
        end

        for tagIdx = 1:numTotalTags
            reset(tgnChannel);
            bxSig{tagIdx} = tgnChannel(bxSig{tagIdx});
        end

        rx = complex(zeros(length(bxSig{1}),1));
        for tagIdx = 1:numTotalTags
            rx = rx+bxSig{tagIdx};
        end

        [rxFromTags,~,~] = func_awgn(rx,snr,'measured');
        ofdmDemod = survey_ConcurScatter_funcReceiver(...
            rxFromTags(ind.HTData(1):ind.HTData(2)),cfgHT,1);
        receivedTagSymbols = ofdmDemod(:,2:1+numTagData);

        [~,originalSymbols] = ...
            survey_ConcurScatter_funcOFDMSymDerived(txPSDU,cfgHT);
        excitationSymbols = originalSymbols(:,2:1+numTagData);

        % The decoder reserves the expanded reference region. The reference
        % estimates of added tags therefore include interference from the
        % initial tags' early payload symbols.
        [decodedSymbols,H_est] = ... %#ok<NASGU>
            survey_ConcurScatter_funcChannelAwareDecode(...
            receivedTagSymbols,excitationSymbols,numTotalTags,...
            seqLenForEstChannel,num_circshift1,num_circshift2);
        decodedBits = pskdemod(decodedSymbols,tag_modulation);

        initialOriginalBits = ...
            pskdemod(initialPayloadTagData,tag_modulation).';
        for tagIdx = 1:numTags
            % Symbols sent by an initial tag while added tags transmit their
            % references are outside the decoder's common payload region.
            % As in Mecha, model these unavailable bits as random guesses.
            unavailableBits = randi(...
                [0,1],numInitialPayloadSymbolsLost,1);
            recoveredInitialBits = ...
                [unavailableBits;decodedBits(:,tagIdx)];
            numBitErrs(i,tagIdx) = numBitErrs(i,tagIdx)+...
                biterr(initialOriginalBits(:,tagIdx),...
                recoveredInitialBits);
        end

        if numAddTags>0
            addedOriginalBits = ...
                pskdemod(addedPayloadTagData,tag_modulation).';
            for tagIdx = numTags+1:numTotalTags
                addedTagIdx = tagIdx-numTags;
                numBitErrs(i,tagIdx) = numBitErrs(i,tagIdx)+...
                    biterr(addedOriginalBits(:,addedTagIdx),...
                    decodedBits(:,tagIdx));
            end
        end

        n = n+1;
    end

    for tagIdx = 1:numTotalTags
        % Use the initial frame's payload length for every tag, matching the
        % normalization in Mecha's dynamic-tag-variation experiment.
        berEst(i,tagIdx) = numBitErrs(i,tagIdx)/...
            (numInitialPayload*maxNumPackets);
    end
end

aaa = 1; %#ok<NASGU>



