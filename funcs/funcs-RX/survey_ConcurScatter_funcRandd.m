function out = survey_ConcurScatter_funcRandd(len,modd)

out = round(rand(1,len)*(modd-1));
out = pskmod(out,modd);

end

