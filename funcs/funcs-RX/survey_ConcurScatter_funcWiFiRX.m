function [rxPSDU,flag] = survey_ConcurScatter_funcWiFiRX(rx,cfgHT)

fs = wlanSampleRate(cfgHT); % Get the baseband sampling rate
ind = wlanFieldIndices(cfgHT); % Indices for accessing each field within the time-domain packet


flag = 1; % 1 -> no packet loss; 0 -> packet loss

coarsePktOffset = wlanPacketDetect(rx,cfgHT.ChannelBandwidth,0,1e-6); % Packet detect and determine coarse packet offset

lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:); % Extract L-STF and perform coarse frequency offset correction
coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgHT.ChannelBandwidth);
rx = helperFrequencyOffset(rx,fs,-coarseFreqOff);

% Extract the non-HT fields and determine fine packet offset
nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
    cfgHT.ChannelBandwidth);

% Determine final packet offset
pktOffset = coarsePktOffset+finePktOffset;
        
if pktOffset>15
    flag = 0;
    rxPSDU = 0;
    return;
end
        
% Extract L-LTF and perform fine frequency offset correction
lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cfgHT.ChannelBandwidth);
rx = helperFrequencyOffset(rx,fs,-fineFreqOff);

% Extract HT-LTF samples from the waveform, demodulate and perform
% channel estimation
htltf = rx(pktOffset+(ind.HTLTF(1):ind.HTLTF(2)),:);
htltfDemod = wlanHTLTFDemodulate(htltf,cfgHT);
chanEst = wlanHTLTFChannelEstimate(htltfDemod,cfgHT);

% Extract HT Data samples from the waveform
htdata = rx(pktOffset+(ind.HTData(1):ind.HTData(2)),:);

% Estimate the noise power in HT data field
nVarHT = htNoiseEstimate(htdata,chanEst,cfgHT);

% Recover the transmitted PSDU in HT Data
rxPSDU = wlanHTDataRecover(htdata,chanEst,nVarHT,cfgHT);

end

