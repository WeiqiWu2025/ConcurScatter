function [decodedSymbols,H_est] = ...
    survey_ConcurScatter_funcChannelAwareDecode(...
    receivedSymbols,excitationSymbols,numTags,seqLenForEstChannel,...
    numCircshiftInteger,numCircshiftFractional)
%survey_ConcurScatter_funcChannelAwareDecode Estimate channels and decode.
%   RECEIVEDSYMBOLS and EXCITATIONSYMBOLS use the full received FFT-bin
%   order. H_EST has dimensions numSubcarriers-by-numTags, where
%   H_EST(n,i) is tag i's estimated channel on received FFT bin n.
%
%   The cyclic shift is included in each reconstructed pattern P_i. H_EST
%   remains indexed by physical received FFT bins and is not shifted again.

numSubcarriers = size(receivedSymbols,1);
numTagData = size(receivedSymbols,2);
numReferenceSymbols = seqLenForEstChannel*numTags;
numPayload = numTagData-numReferenceSymbols;

if size(excitationSymbols,1)~=numSubcarriers || ...
        size(excitationSymbols,2)~=numTagData
    error('ConcurScatter:SymbolDimensionMismatch',...
        'Received and excitation OFDM-symbol matrices must have equal size.');
end
if numPayload<1
    error('ConcurScatter:MissingPayload',...
        'At least one payload OFDM symbol is required.');
end

% Mecha-style per-tag, per-subcarrier complex LS estimation. During tag
% i's reference field only tag i is active, but it retains its assigned
% ConcurScatter cyclic shift.
H_est = complex(zeros(numSubcarriers,numTags));
for tagIdx = 1:numTags
    refCols = (tagIdx-1)*seqLenForEstChannel+...
        (1:seqLenForEstChannel);
    receivedReference = receivedSymbols(:,refCols);
    excitationReference = excitationSymbols(:,refCols);

    referencePattern = survey_ConcurScatter_funcReconstructPattern(...
        excitationReference,numCircshiftInteger(tagIdx),...
        numCircshiftFractional(tagIdx));

    numerator = sum(conj(referencePattern).*receivedReference,2);
    denominator = sum(abs(referencePattern).^2,2);
    valid = denominator>1e-12;
    H_est(valid,tagIdx) = numerator(valid)./denominator(valid);
end

payloadCols = numReferenceSymbols+1:numTagData;
receivedPayload = receivedSymbols(:,payloadCols);
excitationPayload = excitationSymbols(:,payloadCols);

correlationVector = complex(zeros(numTags,numPayload));
channelAwarePatterns = ...
    complex(zeros(numTags,numSubcarriers,numPayload));

for tagIdx = 1:numTags
    idealPattern = survey_ConcurScatter_funcReconstructPattern(...
        excitationPayload,numCircshiftInteger(tagIdx),...
        numCircshiftFractional(tagIdx));
    channelAwarePattern = H_est(:,tagIdx).*idealPattern;

    channelAwarePatterns(tagIdx,:,:) = channelAwarePattern;
    correlationVector(tagIdx,:) = sum(...
        receivedPayload.*conj(channelAwarePattern),1);
end

correlationMatrix = complex(zeros(numTags,numTags,numPayload));
for templateIdx = 1:numTags
    for tagIdx = 1:numTags
        correlationMatrix(templateIdx,tagIdx,:) = sum(...
            channelAwarePatterns(tagIdx,:,:).*...
            conj(channelAwarePatterns(templateIdx,:,:)),2);
    end
end

decodedSymbols = complex(zeros(numPayload,numTags));
for payloadIdx = 1:numPayload
    currentMatrix = reshape(...
        correlationMatrix(:,:,payloadIdx),numTags,numTags);
    currentVector = correlationVector(:,payloadIdx);
    decodedSymbols(payloadIdx,:) = ...
        (pinv(currentMatrix)*currentVector).';
end
end
