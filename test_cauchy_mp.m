rng(0)
addpath('include', 'external');

digs = [32 64 128 256];
n = 100;

delta = 10;
X = gallery('cauchy',n);

mmax = 200;
mp.Digits(32);
Xpowernorm = zeros(1,mmax,'mp');
for i=1:mmax
    Xpowernorm(i) = norm(mp(X)^i,1);
end
num_digs = length(digs);
data = zeros(num_digs, 5);

fileID = fopen('../figs/data_Cauchymat.txt','w');
for i=1:num_digs 
    % fprintf('\n* %d digits...\n', digs(i));
    u = 10^(-digs(i));
    ii = 1:mmax;
    ii = mp(ii);
    tt = Xpowernorm ./ factorial(ii);
    m = sum(tt>u) + 1; % determine the largest tested degree
    [p, s, Cp, di] = PS_mp(X, m, u, delta);
    P = PS_ref(X, m, digs(i)*2);
    err = double(norm(P-p,1)/norm(P,1));
    r = floor(m/s);
    data(i,:) = [err m s r Cp];
    fprintf(fileID,'\n%12s %6s %6s %6s %6s\n','error','m', 's', 'r', 'Cp');
    fprintf(fileID,'%12.2e %6d %6d %6d %6.2f\n', data(i,:));
    fprintf(fileID,'decimal digits:');
    for j=1:r
        fprintf(fileID,'%5d', di(j));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);