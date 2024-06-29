addpath('external','include',"testmats");

format compact
warning off

rng(0);

sizes = [200 500 1000];
nsizes = length(sizes);

deci_digits = [64 256];
num_digits = length(deci_digits);

[~, num_mats] = expm_testmats_size_n(); % number of matrices of variable size

for i=1:nsizes
    n = sizes(i);
    for j=1:num_digits
        digits = deci_digits(j);
        % use mp class to define u, as it can be flushed to zero in double precision
        mp.Digits(34);
        u = 10^(mp(-digits)); 
        fixps_time = zeros(num_mats,1);  % fixed-precision PS.
        mixps_time = zeros(num_mats,4);  % mixed-precision PS.
        sparse_density = zeros(num_mats,1);
        fprintf('\n* Matrix size %d, %d digits \n', n, digits);
        for k=1:1:num_mats
            A = expm_testmats_size_n(k, n);
            sparse_density(k) = nnz(A) / numel(A);
            time = tic;
            fprintf('\n* Matrix id: %d, get params...', k);
            [scal, m] = expm_params_taylor_ap(double(A), 'precision', digits, 'abserr', false);
            fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
            X = A/2^scal;

            time = tic;
            fprintf('Running mixed-precision PS... ');
            [~, ~, ~, ~ , mixps_time(k,:)] = mpps_exp_taylor_ap(X, m, u);
            fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
            
            time = tic;
            fprintf('Running fixed-precision PS... ');
            [~, fixps_time(k)] = fixps_exp_taylor_ap(X, m, digits);
            fprintf('done in [%0.2f minutes]. \n', toc(time)/60);
        end
        time_ratio = sum(mixps_time, 2) ./ fixps_time;
        [~, perm] = sort(sparse_density);
        % save the data for different precisions
        dataname = sprintf('data/exp_taylor_ap_time_%d_%04d.mat', n, digits);
        save(dataname, 'n', 'digits', 'fixps_time', 'mixps_time', 'perm', 'sparse_density', 'time_ratio', 'num_mats');
    end
end

%% load the data and plot

for j=1:num_digits
    digits = deci_digits(j);
    figure
    plot(1:num_mats, ones(1,num_mats), '--', 'LineWidth', 2);
    hold on
    for i=1:nsizes
        n = sizes(i);
        dataname = sprintf('data/exp_taylor_ap_time_%d_%04d.mat', n, digits);
        load(dataname);
        plot(1:num_mats, time_ratio(perm), '-o', 1:num_mats, sparse_density(perm), '-v',...
            'LineWidth',1.5 ,'MarkerSize',3)
    end
    hold off
    timelegend1 = sprintf('time[n=%d]', sizes(1));
    timelegend2 = sprintf('time[n=%d]', sizes(2));
    timelegend3 = sprintf('time[n=%d]', sizes(3));
    sparselegend1 = sprintf('spar[n=%d]', sizes(1));
    sparselegend2 = sprintf('spar[n=%d]', sizes(2));
    sparselegend3 = sprintf('spar[n=%d]', sizes(3));
    legend('ratio=1', timelegend1, timelegend2, timelegend3,...
        sparselegend1, sparselegend2, sparselegend3,...
        'interpreter', 'latex', 'Location', 'NE', 'FontSize', 16);
    set(gca,'linewidth',1.2)
    set(gca,'fontsize',16)
    xlim([0,75]);
    xticks(0:15:75);
    figname = sprintf('data/exp_taylor_ap_time_%04d.eps', digits);
    exportgraphics(gca, figname, 'ContentType', 'vector');
end
