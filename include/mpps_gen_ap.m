function [p, s, cr] = mpps_gen_ap(X, m, coef, u)
%MPPS_GEN_AP    The mixed-precision paterson-stockmeyer algorithm for 
% general polynomials of matrices in arbitrary precision.

% The input COEF is a (m+1)-vector which stores the coefficients of the
% matrix polynomial up to degree m, and u is the unit roundoff of the 
% working precision.
% P is the computed polynomial, the parameter in the PS scheme
% is S = ceil(sqrt(m)), and CR is the complexity reduction.


delta = 10; % precision switching parameter, delta = 10 is recommended

n = size(X,1); % size of the matrix
s = ceil(sqrt(m));  % the (default) parameter s in the PS scheme


mp.Digits(34);
u = mp(u);
d_u = double(round(log10(1/u),0));
mp.Digits(d_u);

Xsq_powers = zeros(n,n,s,'mp');
Xsq_powers(:,:,1) = eye(n,'mp');
Xsq_powers(:,:,2) = X;

[B0, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u);
normB0 = normest1(double(B0));
normY = normest1(double(Y));

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

% compute B_i using Xsq_powers, estimate ||B_i|| and decide whether to use
% lower precisions
for j=1:r
    Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);

    mp.Digits(34);
    Bsq_coefmat_norm(j+1) = normest1(mp(Bsq_coefmat(:,:,j+1))); 
    usq_prec(j+1) = max(u, Bsq_coefmat_norm(1)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j))*u);
    mp.Digits(d_u);
    if usq_prec(j+1)>=delta*u
        nu = j;
        break;
    end
end

for j=nu+1:r
    Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
    mp.Digits(34);
    Bsq_coefmat_norm(j+1) = normest1(mp(Bsq_coefmat(:,:,j+1))); 
    usq_prec(j+1) = Bsq_coefmat_norm(nu)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j-nu+1))*u;
    mp.Digits(d_u);
end
usq_prec = min(usq_prec, 0.1); % use at least 1 digit in the computation

% compute the polynomial p using Horner's method
p = Bsq_coefmat(:,:,r+1);
for k=r:-1:1
    d_ui = double(round(log10(1/usq_prec(k+1)), 0));
    mp.Digits(d_ui);
    p = mp(p)*mp(Y);
    d_ui = double(round(log10(1/usq_prec(k)), 0));
    mp.Digits(d_ui);
    p = mp(p) + mp(Bsq_coefmat(:,:,k));
end

divec = double( round(log10(  1./ usq_prec(1:(r+1))  ), 0) );
cr = 1 - ((s-1)*divec(1) + sum(divec)-divec(1)) / ((s+r-1)*divec(1)); % the complexity reduction


%% -----------------------       SUBFUNCTIONS       -----------------------
function [q, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u)
% evaluate the coefficient matrix B_0(X) in precision u

mp.Digits(d_u);

q = coef(1)*Xsq_powers(:,:,1);

for i=1:s-1
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
    q = q + coef(i+1)*Xsq_powers(:,:,i+1);
end
Y = Xsq_powers(:,:,s)*Xsq_powers(:,:,2);
end

function q = eval_coef_mat(Xsq_powers, s, d_u, coef_mat_index)
% form B_i using the powers from Xsq_powers

mp.Digits(d_u);

if coef_mat_index<r
    % form B_i(X), i=1:r-1
    num_term = s-1;
    coef_B_i = coef((coef_mat_index*s+1):((coef_mat_index+1)*s))';
else
    % form B_r(X)
    num_term = m-s*r;
    coef_B_i = coef((coef_mat_index*s+1):(coef_mat_index*s+num_term+1))';
end
reshaped_vector = reshape(coef_B_i,[1,1,num_term+1]);
q = sum(bsxfun(@times, Xsq_powers(:,:,1:(num_term+1)),reshaped_vector),3);
end


end