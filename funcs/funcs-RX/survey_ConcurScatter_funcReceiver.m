function ofdmDemod = survey_ConcurScatter_funcReceiver(rxHTData,cfgHT,numSTS)


% Get OFDM related parameters
[cfgOFDM,dataInd,pilotInd] = wlan.internal.wlanGetOFDMConfig(cfgHT.ChannelBandwidth, cfgHT.GuardInterval, 'HT', numSTS);
ofdmSymOffset = 0.75;

FFTLen    = cfgOFDM.FFTLength;
CPLen     = cfgOFDM.CyclicPrefixLength;
numRx     = size(rxHTData, 2); 
symOffset = round(ofdmSymOffset * CPLen);

% Remove cyclic prefix
numSym = size(rxHTData,1)/(FFTLen + CPLen);
inputIn3D = reshape(rxHTData, [(FFTLen + CPLen) numSym numRx]);
postCPRemoval = inputIn3D([CPLen+1:FFTLen+symOffset, symOffset+1:CPLen], :, :);

% Denormalization
postCPRemoval = postCPRemoval / cfgOFDM.NormalizationFactor;

% FFT
postFFT = fft(postCPRemoval, [], 1);

% FFT shift
if isreal(postFFT)
    postShift = complex(fftshift(postFFT, 1), 0);
else
    postShift = fftshift(postFFT,1);
end

% Phase rotation on frequency subcarriers
postShift = bsxfun(@rdivide, postShift, cfgOFDM.CarrierRotations);

ofdmDemod = postShift;

end

