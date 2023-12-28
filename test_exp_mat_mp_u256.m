addpath('include', 'external');
rng(0);


n = 100;
num_digs = 256;


delta = 10;

u = 10^(-num_digs);
[~, num_total_mats] = mygallery();

num_mats = 97;

rr = zeros(num_mats,1);
ee = zeros(num_mats,1);
ee_fixprec = zeros(num_mats,1);
n_mats = zeros(num_mats,1);
complex_reduc = zeros(num_mats,1);


ii = 1;
i = 1;
while i<=num_total_mats
    % fprintf('\n* Matrix id: %d\n', i);
    A = mygallery(i,n);
    while ishermitian(A) % exclude Hermitian matrices
        i = i + 1;
        A = mygallery(i,n);
    end
    [~, scal, m] = expm_mp(double(A),...        
                                                'precision',num_digs,...
                                                'algorithm', 'transfree',...
                                                'abserr', false,...
                                                'approx', 'taylor');
    X = A/2^scal;
    n_mats(ii) = size(X,1);
    [p, s, complex_reduc(ii)] = PS_mp(X, m, u, delta);
    p_fixprec = PS_fix_mp(X, m, num_digs);
    P = PS_ref(X, m, 400);
    rr(ii) = floor(m/s);
    normP = norm(P,1);
    ee(ii) = double(norm(P-p,1)/normP);
    ee_fixprec(ii) = double(norm(P-p_fixprec,1)/normP);
    i = i + 1;
    ii = ii + 1;
end
[~, perm] = sort(n_mats.*rr);

ee_fixprec_scal = scal_small_error(ee_fixprec, num_digs);
ee_scal = scal_small_error(ee, num_digs);


% save the data
dataname = sprintf('data/exp_mat_mp_u256.mat');
save(dataname, 'n_mats', 'rr', 'perm', 'u', 'ee_scal', 'ee_fixprec_scal', 'complex_reduc', 'num_mats');
   
% load the data
% dataname = sprintf('data/exp_mat_mp_u256.mat');
% load(dataname)

figure

semilogy(1:num_mats, rr(perm).*n_mats(perm)*u,'-', 1:num_mats, ee_scal(perm),'v', 1:num_mats, ee_fixprec_scal(perm),'o', ...
    1:num_mats, u*ones(num_mats,1),'--','LineWidth',2,'MarkerSize',8)

mycolors = [0 0 0; 0.2300 0.4800 0.3400; 1 0.4900 0; 0.3010 0.7450 0.9330];
ax = gca; 
ax.ColorOrder = mycolors;


legend('$rnu$', '$\epsilon_{v}$', '$\epsilon_{f}$', '$u$', 'interpreter', 'latex', 'Location', 'NE', 'FontSize', 16);
set(gca,'linewidth',1.2)
set(gca,'fontsize',16)


xlim([0,num_mats+1]);
xticks([0:15:90 num_mats]);

ylim([1e-257 1e-247])
yticks(10.^(-257:2:-247))

exportgraphics(gca, '../figs/err-mat_n100_u256.eps', 'ContentType', 'vector');


figure
bar(complex_reduc(perm))
set(gca,'linewidth',1.2)
set(gca,'fontsize',16)
xlim([0,num_mats+1]);
xticks([0:15:90 num_mats]);
ynum=[cellstr(num2str(get(gca,'ytick')'*100))];
pct = char(ones(size(ynum,1),1)*'%'); % Create a vector of '%' signs.
new_yticks = [char(ynum),pct]; % Append the '%' signs after the percentage values.
yticklabels(new_yticks);

exportgraphics(gca, '../figs/cmplxreduc-mat_n100_u256.eps', 'ContentType', 'vector');
