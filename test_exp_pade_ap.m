addpath('external','include',"testmats");

format compact
warning off

rng(0);

n = 100;
num_digs = 64;
u = 10^(-num_digs);

yaxis_lim =  [1e-65 1e-59];
yaxis_ticks = 10.^(-65:1:-59);

padetype = ["numerator", "denominator"];

[~, num_mats] = expm_testmats(); % total number of matrices 

for i=1:length(padetype)
    padepoly = padetype(i); % numerator or denominator polynomial in Pade
    rr = zeros(num_mats,1);
    err_mpps = zeros(num_mats,1);
    err_fixps = zeros(num_mats,1);
    size_mats = zeros(num_mats,1);
    complex_reduc = zeros(num_mats,1);
    fprintf('\n* Tesing %s polyn in Pade, Matrix size %d \n', padepoly, n);
    for k=1:num_mats
        A = expm_testmats(k, n);
        time = tic;
        fprintf('\n* Matrix id: %d, get Pade params...', k);
        [scal, m] = expm_params_pade_ap(double(A), 'precision', num_digs, 'abserr', false);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);

        X = A/2^scal;
        size_mats(k) = size(X, 1);

        time = tic;
        fprintf('Computing Pade coefficients... ');
        coef = comput_coef_pade_ap(m, num_digs, padepoly);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
        
        time = tic;
        fprintf('Running mixed-precision PS... ');
        [p_mixprec, s, complex_reduc(k)] = mpps_gen_ap(X, m, coef, u);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);

        time = tic;
        fprintf('Running fixed-precision PS... ');
        p_fixprec = fixps_gen_ap(X, m, coef, num_digs);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
    
        time = tic;
        fprintf('Calculating reference PS... ');
        p_ref = refps_gen(X, m, coef, num_digs);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
        
        rr(k) = floor(m/s);
        norm_p_ref = norm(p_ref,1);
        err_mpps(k) = double(norm(p_ref-p_mixprec,1)/norm_p_ref);
        err_fixps(k) = double(norm(p_ref-p_fixprec,1)/norm_p_ref);
    end
    [~, perm] = sort(size_mats.*rr);
    err_fixps_scal = scal_small_error(err_fixps, num_digs);
    err_mpps_scal = scal_small_error(err_mpps, num_digs);
    % save the data for different types of Pade polynomials
    dataname = sprintf('data/exp_pade_ap_error_%d_%04d_%s.mat', n, num_digs, padepoly);
    save(dataname, 'n', 'num_digs', 'err_mpps_scal', 'err_fixps_scal', ...
            'complex_reduc', 'perm', 'rr', 'num_mats','u');
end

%% load the data and plot

for i=1:length(padetype)
    padepoly = padetype(i);
    figure
    dataname = sprintf('data/exp_pade_ap_error_%d_%04d_%s.mat', n, num_digs, padepoly);
    load(dataname);
    semilogy(1:num_mats, rr(perm).*size_mats(perm)*u, '-',...
        1:num_mats, err_mpps_scal(perm), 'v',...
        1:num_mats, err_fixps_scal(perm),'o', ...
        1:num_mats, u*ones(num_mats,1),  '--', 'LineWidth', 2,'MarkerSize', 8)
    mycolors = [0 0 0; 0.2300 0.4800 0.3400; 1 0.4900 0; 0.3010 0.7450 0.9330];
    ax = gca; 
    ax.ColorOrder = mycolors;
    legend('$rnu$', '$\epsilon_{v}$', '$\epsilon_{f}$', '$u$', ...
        'interpreter', 'latex', 'Location', 'NE', 'FontSize', 16);
    set(gca,'linewidth',1.2)
    set(gca,'fontsize',16)
    xlim([0,num_mats+1]);
    xticks([0:20:120 num_mats]);
    ylim(yaxis_lim)
    yticks(yaxis_ticks)
    figname1 = sprintf('data/exp_pade_ap_error_%d_%04d_%s.eps', n, num_digs, padepoly);
    exportgraphics(gca, figname1, 'ContentType', 'vector');
    figure
    bar(complex_reduc(perm))
    set(gca,'linewidth',1.2)
    set(gca,'fontsize',16)
    xlim([0,num_mats+1]);
    xticks([0:20:120 num_mats]);
    ynum=[cellstr(num2str(get(gca,'ytick')'*100))];
    pct = char(ones(size(ynum,1),1)*'%'); % Create a vector of '%' signs.
    new_yticks = [char(ynum),pct]; % Append the '%' signs after the percentage values.
    yticklabels(new_yticks);
    figname2 = sprintf('data/exp_pade_ap_cmplxreduc_%d_%04d_%s.eps', n, num_digs, padepoly);
    exportgraphics(gca, figname2, 'ContentType', 'vector');
end

%% subfunction

function coef = comput_coef_pade_ap(m, num_digs, options)
% compute the vector of Pade coefficients for the matrix exponential in
% arbitrary precision.

mp.Digits(num_digs);
m = mp(m);
l_num = mp(0:1:m);
if isequal(options,'numerator')
    coef = (factorial(m)/factorial(2 * m)) ./...
                (factorial(m-l_num).*factorial(l_num)) .*...
                factorial(2*m - l_num);
end
if isequal(options,'denominator')
    coef = (-1).^l_num .* ...
        (factorial(m)/factorial(2 * m)) ./ ...
                (factorial(m-l_num).*factorial(l_num)) .* ...
                factorial(2*m - l_num);
end
end