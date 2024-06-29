function [p, s, cr, divec, time] = mpps_exp_taylor_ap(X, m, u)
%PS_EXP_TAYLOR_AP   The mixed-precision paterson-stockmeyer algorithm for 
% the Taylor approximants of the matrix exponential in arbitrary precision. 
% 
% The degree of the Taylor approximant is specified by M, 
% U is the unit roundoff of the working precision,
% The algorithm can also report the computational time TIME. P is the 
% computed polynomial, the parameter in the PS scheme is S = ceil(sqrt(m)).
% CR is the complexity reduction, and DIVEC returns the vector of used 
% decimal digits in the mixed-precision computation.



% time(1): Computing the required matrix powers in precision u.
% time(2): Total norm estimation time.
% time(3): The Horner's stage: matmul and matadd in mixed precisions.
% time(4): Total time for assembling the coefficient matrices.
time = zeros(1,4);

delta = 10; % precision switching parameter, delta = 10 is recommended

n = size(X,1); % size of the matrix
s = ceil(sqrt(m));  % the (default) parameter s in the PS scheme

mp.Digits(34);
u = mp(u);
d_u = double(round(log10(1/u),0));

mp.Digits(d_u);

Xsq_powers = zeros(n,n,s,'mp');
Xsq_powers(:,:,1) = eye(n);
Xsq_powers(:,:,2) = X;

[B0, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u);

normest_time = tic;
normB0 = normest1(double(B0));
normY = normest1(double(Y));
time(2) = toc(normest_time);

r = floor(m/s);
nu = r+1;
Bsq_coefmat = zeros(n,n,r+1,'mp');
Bsq_coefmat(:,:,1) = B0;

% use quadruple precision for norm estimation to avoid over- and underflow
mp.Digits(34);
Bsq_coefmat_norm = zeros(1,r+1,'mp');
usq_prec = zeros(1,r+1,'mp');
mp.Digits(d_u);

Bsq_coefmat_norm(1) = normB0;
usq_prec(1) = u;
for j=1:r
    assemb_time = tic;
    Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
    time(4) = time(4) + toc(assemb_time);
    mp.Digits(34);

    normest_time = tic;
    Bsq_coefmat_norm(j+1) = normest1(mp(Bsq_coefmat(:,:,j+1))); 
    time(2) = time(2) + toc(normest_time);

    usq_prec(j+1) = max(u, Bsq_coefmat_norm(1)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j))*u);
    mp.Digits(d_u);
    if usq_prec(j+1)>=delta*u
        nu = j;
        break;
    end
end
for j=nu+1:r

    assemb_time = tic;
    Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
    time(4) = time(4) + toc(assemb_time); 

    mp.Digits(34);

    normest_time = tic;
    Bsq_coefmat_norm(j+1) = normest1(mp(Bsq_coefmat(:,:,j+1))); 
    time(2) = time(2) + toc(normest_time);

    usq_prec(j+1) = Bsq_coefmat_norm(nu)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j-nu+1))*u;
    mp.Digits(d_u);
end
usq_prec = min(usq_prec, 0.1); % use at least 1 digit in the computation

% compute the polynomial p using Horner's method
horner_time = tic;
p = Bsq_coefmat(:,:,r+1);
for k=r:-1:1
    d_ui = double(round(log10(1/usq_prec(k+1)), 0));
    mp.Digits(d_ui);
    p = mp(p)*mp(Y);
    d_ui = double( round(log10(1/usq_prec(k)), 0));
    mp.Digits(d_ui);
    p = mp(p) + mp(Bsq_coefmat(:,:,k));
end
time(3) = toc(horner_time);

divec = double( round(log10(  1./ usq_prec(1:(r+1))  ), 0) );
cr = 1 - ((s-1)*divec(1) + sum(divec)-divec(1)) / ((s+r-1)*divec(1)); % the complexity reduction


%% -----------------------       SUBFUNCTIONS       -----------------------
function [q, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u)
% evaluate the coefficient matrix B_0(X) in precision u

mp.Digits(d_u);

matpow_time = tic;
for i=2:s-1
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*X;
end
Y = Xsq_powers(:,:,s)*X;
time(1) = toc(matpow_time);

assemb_time = tic;
q = Xsq_powers(:,:,s) / (s-1);
for i = s-2:-1:1
    q = (q + Xsq_powers(:,:,i+1)) / i;
end
q(1:n+1:n^2) = q(1:n+1:n^2) + 1;
time(4) = time(4) + toc(assemb_time);
end

function q = eval_coef_mat(Xsq_powers, s, d_u, coef_mat_index)
% form B_i using the powers of X from Xsq_powers

mp.Digits(d_u);
if coef_mat_index<r
    % form B_i(X), i=1:r-1
    q = Xsq_powers(:,:,s) / (s*coef_mat_index + s-1);
    for i = s-2:-1:0
        q = (q + Xsq_powers(:,:,i+1)) / (s*coef_mat_index + i);
    end
    q = q / factorial(mp(s*coef_mat_index-1));
else
    % form B_r(X), num_term = m-s*r;
    q = Xsq_powers(:,:,m-s*r+1) / m;
    for i = m-s*r-1:-1:0
        q = (q + Xsq_powers(:,:,i+1)) / (s*r + i);
    end
    q = q / factorial(mp(s*r-1));
end
end

end