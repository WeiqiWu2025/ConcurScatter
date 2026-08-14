function est_H = survey_ConcurScatter_funcChannelEstimator(ofdm_demd,ofdm_reference)

tmpLL = size(ofdm_demd,1);
for idx1 = 1:tmpLL
    tmp_H_real = funcLSEstimator(ofdm_reference(idx1,:)',real(ofdm_demd(idx1,:))');
    tmp_H_imag = funcLSEstimator(ofdm_reference(idx1,:)',imag(ofdm_demd(idx1,:))');
    tmp_H = tmp_H_real + 1i*tmp_H_imag;
    est_H(idx1) = tmp_H;
end

end

