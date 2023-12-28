rng(0)
addpath('include', 'external');


e = exp(1);
n = 50;

delta = 10;

num_digs = 64;
u = 10^(-num_digs);

mp.Digits(100);
ii = 1:300;
ii = mp(ii);
tt = ceil(sqrt(ii)).^ii ./ factorial(ii) ./ exp(ii);

mmax = sum(tt>u) + 1; % determine the largest tested degree



mmin = 5;
mgap = 3;
mm = mmin:mgap:mmax;

num = length(mm);


ss = ceil(sqrt(mm));
rr = zeros(num,1);
sigma = zeros(num,1);
for i=1:num
    sigma(i) = ss(i)/e;
end


% now X = rand(n)
X = rand(n);

ee = zeros(num,1);
ee_fixprec = zeros(num,1);

for i=1:num
    % fprintf('\n* X = rand(n), m = %d\n', mm(i));
    Xnorm = norm(X,1);
    X = X/Xnorm*sigma(i);
    [p, s] = PS_mp(X, mm(i), u, delta);
    p_fixprec = PS_fix_mp(X, mm(i), num_digs);
    P = PS_ref(X, mm(i));
    normP = norm(P,1);
    rr(i) = floor(mm(i)/s);
    ee(i) = double(norm(P-p,1)/normP);
    ee_fixprec(i) = double(norm(P-p_fixprec,1)/normP);
end

ee_fixprec_scal = scal_small_error(ee_fixprec, num_digs);
ee_scal = scal_small_error(ee, num_digs);

% save the data
dataname = sprintf('data/exp_m_mp_rand.mat');
save(dataname, 'mm', 'rr', 'n', 'u', 'ee_scal', 'ee_fixprec_scal', 'num');
   
% load the data
% dataname = sprintf('data/exp_m_mp_rand.mat');
% load(dataname)


figure
semilogy(mm, rr*n*u, '-', ...
    mm, ee_scal,'-v', ...
    mm, ee_fixprec_scal,'-o', ...
    mm, u*ones(num,1), '--', ...
    'LineWidth',2,'MarkerSize',8)
mycolors = [0 0 0; 0.2300 0.4800 0.3400; 1 0.4900 0; 0.3010 0.7450 0.9330];
ax = gca; 
ax.ColorOrder = mycolors;


legend('$rnu$', '$\epsilon_v$', '$\epsilon_f$', '$u$', 'interpreter', 'latex', 'Location', 'NW', 'FontSize', 16);
set(gca,'linewidth',1.2)
set(gca,'fontsize',16)
xlabel('$m$','interpreter','latex','FontWeight','normal','fontsize',16)
xlim([mmin mmax])
xticks(mmin:2*mgap:mmax)

ylim([1e-65 1e-60])
yticks(10.^(-65:1:60))

% ylim([5e-17 1e-12])
exportgraphics(gca, '../figs/err-m_rand_n50_u64.eps', 'ContentType', 'vector');



% now X = randn(n)
X = randn(n);


ee = zeros(num,1);
ee_fixprec = zeros(num,1);

for i=1:num
    % fprintf('\n* X = randn(n), m = %d\n', mm(i));
    Xnorm = norm(X,1);
    X = X/Xnorm*sigma(i);
    [p, s] = PS_mp(X, mm(i), u, delta);
    p_fixprec = PS_fix_mp(X, mm(i), num_digs);
    P = PS_ref(X, mm(i));
    normP = norm(P,1);
    rr(i) = floor(mm(i)/s);
    ee(i) = double(norm(P-p,1)/normP);
    ee_fixprec(i) = double(norm(P-p_fixprec,1)/normP);
end

ee_fixprec_scal = scal_small_error(ee_fixprec, num_digs);
ee_scal = scal_small_error(ee, num_digs);

% save the data
dataname = sprintf('data/exp_m_mp_randn.mat');
save(dataname, 'mm', 'rr', 'n', 'u', 'ee_scal', 'ee_fixprec_scal', 'num');
   
% load the data
% dataname = sprintf('data/exp_m_mp_randn.mat');
% load(dataname)


figure
semilogy(mm, rr*n*u, '-', ...
    mm, ee_scal,'-v', ...
    mm, ee_fixprec_scal,'-o', ...
    mm, u*ones(num,1), '--', ...
    'LineWidth',2,'MarkerSize',8)
mycolors = [0 0 0; 0.2300 0.4800 0.3400; 1 0.4900 0; 0.3010 0.7450 0.9330];

ax = gca; 
ax.ColorOrder = mycolors;


legend('$rnu$', '$\epsilon_v$', '$\epsilon_f$', '$u$', 'interpreter', 'latex', 'Location', 'NW', 'FontSize', 16);
set(gca,'linewidth',1.2)
set(gca,'fontsize',16)
xlabel('$m$','interpreter','latex','FontWeight','normal','fontsize',16)
xlim([mmin mmax])
xticks(mmin:2*mgap:mmax)

ylim([1e-65 1e-60])
yticks(10.^(-65:1:60))

exportgraphics(gca, '../figs/err-m_randn_n50_u64.eps', 'ContentType', 'vector');