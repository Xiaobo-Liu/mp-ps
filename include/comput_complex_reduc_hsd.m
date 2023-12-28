function Cp = comput_complex_reduc_hsd(usq_prec, s, m)
% compute complexity reduction in percentage from the vector USQ_PREC which
% stores the used precisions.

uf = double(2^(-11));
us = double(eps('single')/2);
r = length(usq_prec) - 1;

num_half = sum(usq_prec>=uf);
num_sgl = sum(usq_prec>=us) - num_half;
num_dbl = r - num_sgl - num_half;

s_fix = ceil(sqrt(m));
r_fix = floor(m/s_fix);

Cp = 1 - (4*(s-1+num_dbl)+2*num_sgl+num_half) / (4*(s_fix-1+r_fix));


end