function [preamble] = survey_ConcurScatter_funcGeneratePreamble(Lones_zeroes,numTags)

preamble = [];
for i =1:numTags
    tmp = zeros(Lones_zeroes*numTags,1);
    insertData = ones(Lones_zeroes,1);
    tmp((i-1)*length(insertData)+1:i*length(insertData)) = insertData;
    preamble = [preamble,tmp];
end     

end


