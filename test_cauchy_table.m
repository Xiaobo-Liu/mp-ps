addpath('external','include',"testmats");

format compact
warning off

rng(0);

deci_digits = [32 64 128 256];
length_digits = length(deci_digits);

n = 100;
A = gallery('cauchy',n);

filename = sprintf('tabs/table_cauchy_%d.tex', n);
fileid_cauchy = fopen(filename,'w');
fprintf(fileid_cauchy, ['\\begin{tabularx}{\\textwidth}',...
    '{@{\\extracolsep{\\fill}}crrrrr}\n']);
    fprintf(fileid_cauchy, '\\toprule\n');
    fprintf(fileid_cauchy, ['decimal digits & $m$ & ',...
        '$s$ & $r$ & $(d_1,d_2,\\dots,d_r)$  & $C_r$ \\\\\n']);
    fprintf(fileid_cauchy, '\\midrule \n');
for i=1:length_digits
    num_digs = deci_digits(i);
    fprintf('Num of digits: %d... \n', num_digs);
    [scal, m] = expm_params_taylor_ap(double(A), 'precision', num_digs, 'abserr', false);
    X = A/2^scal;

    u = 10^(-num_digs);
    [p_mixprec, s, cr, divec] = mpps_exp_taylor_ap(X, m, u);

    r = floor(m/s);
    p_ref = refps_exp_taylor(X, m, num_digs);
    norm_p_ref = norm(p_ref,1);
    err_mpps = double(norm(p_ref-p_mixprec,1)/norm_p_ref);
    % all required data obtained, now printing the results
    di_string = sprintf('%g, ', divec(2:r+1));
    di_string = di_string(1:end-2); % remove the extra symbols
    fprintf(fileid_cauchy,...
        '%d & %d & %d & %d & $(%s)$ & %2.1f\\%% \\\\\n', num_digs, m, s, r, di_string, cr*100);
end
fprintf(fileid_cauchy, '\\bottomrule\n');
fprintf(fileid_cauchy, '\\end{tabularx}');
fclose(fileid_cauchy);