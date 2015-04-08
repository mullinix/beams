function F = feuniformforce(q_0,len)
F = (q_0/12).*[6*len; len^2; 6*len; -len^2];
end