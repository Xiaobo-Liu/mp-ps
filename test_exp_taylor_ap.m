addpath('external','include',"testmats");

format compact
warning off

rng(0);

n = 100;

deci_digits = [64 256];

yaxis_lim =  [1e-65 1e-59;
              1e-257 1e-251];
yaxis_ticks = [10.^(-65:1:-59);
               10.^(-257:1:-251);];

length_digits = length(deci_digits);
[~, num_mats] = expm_testmats(); % number of matrices of variable size


for i=1:length_digits
    num_digs = deci_digits(i);
    u = 10^(-num_digs);
    rr = zeros(num_mats,1);
    err_mpps = zeros(num_mats,1);
    err_fixps = zeros(num_mats,1);
    size_mats = zeros(num_mats,1);
    complex_reduc = zeros(num_mats,1);
    fprintf('\n* Test of %d decimal digits, Matrix size %d \n', num_digs, n);
    for k=1:num_mats
        A = expm_testmats(k, n);
        time = tic;
        fprintf('\n* Matrix id: %d, get Taylor params...', k);
        [scal, m] = expm_params_taylor_ap(double(A), 'precision', num_digs, 'abserr', false);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);

        X = A/2^scal;
        size_mats(k) = size(X, 1);
        
        time = tic;
        fprintf('Running mixed-precision PS... ');
        [p_mixprec, s, complex_reduc(k)] = mpps_exp_taylor_ap(X, m, u);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);

        time = tic;
        fprintf('Running fixed-precision PS... ');
        p_fixprec = fixps_exp_taylor_ap(X, m, num_digs);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
    
        time = tic;
        fprintf('Calculating reference PS... ');
        p_ref = refps_exp_taylor(X, m, num_digs);
        fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
        
        rr(k) = floor(m/s);
        norm_p_ref = norm(p_ref,1);
        err_mpps(k) = double(norm(p_ref-p_mixprec,1)/norm_p_ref);
        err_fixps(k) = double(norm(p_ref-p_fixprec,1)/norm_p_ref);
    end
    [~, perm] = sort(size_mats.*rr);
    err_fixps_scal = scal_small_error(err_fixps, num_digs);
    err_mpps_scal = scal_small_error(err_mpps, num_digs);
    % save the data for different working precision enviroments
    dataname = sprintf('data/exp_taylor_ap_error_%d_%04d.mat', n, num_digs);
    save(dataname, 'n', 'num_digs', 'err_mpps_scal', 'err_fixps_scal', ...
            'complex_reduc', 'perm', 'rr', 'num_mats','u');
end

%% load the data and plot

for i=1:length_digits
    num_digs = deci_digits(i);
    figure
    dataname = sprintf('data/exp_taylor_ap_error_%d_%04d.mat', n, num_digs);
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
    ylim(yaxis_lim(i,:))
    yticks(yaxis_ticks(i,:))
    figname1 = sprintf('data/exp_taylor_ap_error_%d_%04d.eps', n, num_digs);
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
    figname2 = sprintf('data/exp_taylor_ap_cmplxreduc_%d_%04d.eps', n, num_digs);
    exportgraphics(gca, figname2, 'ContentType', 'vector');
end