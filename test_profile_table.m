warning off
addpath('include', 'external');
rng(0);
format compact

timing = true;

sizes = [20, 50, 100, 200, 500, 1000];
nsizes = length(sizes);

deci_digits = [64 256];
num_digits = length(deci_digits);

filename = sprintf('tabs/table_profile_%04d_%04d.tex', deci_digits(1), deci_digits(2));
fileid_profile = fopen(filename,'w');
fprintf(fileid_profile, ['\\begin{tabularx}{1.095\\textwidth}',...
        '{@{\\extracolsep{\\fill}}rr|rrrrrrr|rrrrrrr}\n']);
    fprintf(fileid_profile, '\\toprule\n');
    fprintf(fileid_profile, '\\multicolumn{2}{c|}{} & \\multicolumn{7}{c|}{ %d decimal digits} & \\multicolumn{7}{c}{ %d decimal digits} \\\\\n',...
        deci_digits(1), deci_digits(2));
    fprintf(fileid_profile, [' & $n$ & ',...
        '$M_{low}$ & $T_{pow}$ & $T_{est}$ & $T_{hon}$ & $T_{coe}$ & $T_{tot}$ & $T_{fix}$ & ',...
        '$M_{low}$ & $T_{pow}$ & $T_{est}$ & $T_{hon}$ & $T_{coe}$ & $T_{tot}$ & $T_{fix}$\\\\\n']);
    fprintf(fileid_profile, '\\midrule\n');

    mats_id = {'\verb|A|', '\verb|B|', '\verb|C|'};

    nmat = 3;
    for j = 1:nmat
        fprintf('Matrix ID: %1d... \n', j);
        for i = 1:nsizes

            n = sizes(i);
            fprintf('Matrix size = %4d... \n', n);
            switch j
                case 1
                    A = gallery('smoke',n); % sparse
                case 2
                    A = 100 * triu(randn(n),1); % strictly upper-triangular
                case 3
                    A = gallery('lotkin', n); % full
            end
            
            Mtot = zeros(num_digits,1); % total number of matrix products
            Mlow = zeros(num_digits,1); % matmul done in u^{1/2} or higher
            Rlow = zeros(num_digits,1); % ratio of matmul done in u^{1/2} or higher
            Ttot = zeros(num_digits,1);
            Tpow = zeros(num_digits,1);
            Test = zeros(num_digits,1);
            Thon = zeros(num_digits,1);
            Tcoe = zeros(num_digits,1);
            Tfix = zeros(num_digits,1);
            for k=1:num_digits
                [scal, m] = expm_params_taylor_mp(double(A), 'precision', deci_digits(k), 'abserr', false);
                X = A/2^scal;
                mp.Digits(34);
                u = 10^(mp(-deci_digits(k)));
                fprintf('Running mixed-precision PS, digits = %d... \n', deci_digits(k));
                [p, s, ~, divec, mixps_time] = PS_exp_Taylor_mp(X, m, u, timing);
                fprintf('Running fixed-precision PS, digits = %d... \n', deci_digits(k));
                [~, Tfix(k)] = PS_fix_mp(X, m, deci_digits(k), timing);
            
                Mtot(k) = s + floor(m/s) - 1; % total number of matrix products
                Mlow(k) = sum(divec<=deci_digits(k)/2); % matmul done in u^{1/2} or higher
                Rlow(k) = Mlow(k)/Mtot(k)*100;
                Ttot(k) = sum(mixps_time)/100;
                Tpow(k) = mixps_time(1);
                Test(k) = mixps_time(2);
                Thon(k) = mixps_time(3);
                Tcoe(k) = mixps_time(4);
            end
            if (i == 1)
                id_field = mats_id{j};
            else
                id_field = '        ';
            end
            fprintf(fileid_profile,...
                '%s & %3d & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %4.1f & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %2.0f\\%% & %4.1f & %4.1f \\\\\n',...
                id_field, sizes(i),...
                Rlow(1), Tpow(1)/Ttot(1), Test(1)/Ttot(1), Thon(1)/Ttot(1), Tcoe(1)/Ttot(1), Ttot(1)*100, Tfix(1),...
                Rlow(2), Tpow(2)/Ttot(2), Test(2)/Ttot(2), Thon(2)/Ttot(2), Tcoe(2)/Ttot(2), Ttot(2)*100, Tfix(2));

        end
        if (i == nsizes && j ~= nmat)
            fprintf(fileid_profile, '\\midrule\n');
        end
    end
fprintf(fileid_profile, '\\bottomrule\n');
fprintf(fileid_profile, '\\end{tabularx}');