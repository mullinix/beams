function F = ferampedforce(q_0,len)
F = (q_0*len/60).*[9; 2*len; 21; -3*len];
end