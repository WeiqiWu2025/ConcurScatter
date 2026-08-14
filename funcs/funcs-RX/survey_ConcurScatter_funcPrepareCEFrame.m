function [tagData,payloadTagData,seqLenUsed,numReferenceSymbols,numPayload] = ...
    survey_ConcurScatter_funcPrepareCEFrame(...
    numTagData,numTags,requestedSeqLen,tagModulation)
%survey_ConcurScatter_funcPrepareCEFrame Build reference and payload fields.
%   Rows of TAGDATA correspond to tags and columns correspond to OFDM
%   symbols. Each tag receives a disjoint reference field. In tag i's
%   reference field, tag i sends known symbol +1 and all other tags are
%   silent. The remaining symbols carry concurrent PSK payload data.

maxFeasibleSeqLen = floor((numTagData-1)/numTags);
seqLenUsed = min(requestedSeqLen,maxFeasibleSeqLen);

if seqLenUsed<1
    error('ConcurScatter:InsufficientReferenceSymbols',...
        ['The HT-Data field must contain at least one reference symbol ',...
        'per tag and one payload symbol.']);
end

numReferenceSymbols = seqLenUsed*numTags;
numPayload = numTagData-numReferenceSymbols;

tagData = complex(zeros(numTags,numTagData));
for tagIdx = 1:numTags
    refCols = (tagIdx-1)*seqLenUsed+(1:seqLenUsed);
    tagData(tagIdx,refCols) = 1;
end

payloadTagData = reshape(...
    survey_ConcurScatter_funcRandd(...
    numTags*numPayload,tagModulation),numPayload,[]).';
tagData(:,numReferenceSymbols+1:end) = payloadTagData;
end
