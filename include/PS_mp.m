function [p, s, Cp, di_vec] = PS_mp(X, m, u, delta)
% The mixed-precision paterson-stockmeyer algorithm that can use arbitrary
% precision. 
% 
% The degree of the Taylor approximant is specified by m, 
% u is the unit roundoff of the working precision,
% and delta = 10 is recommended.

n = size(X,1); % size of the matrix
s = ceil(sqrt(m));  % the (default) parameter s in the PS scheme
sfactorial = factorial(s);
e = exp(1);

mp.Digits(34);
d_u = double(round(log10(1/mp(u)),0));
mp.Digits(d_u);

Xsq_powers = zeros(n,n,s,'mp');
Xsq_powers(:,:,1) = eye(n,'mp');
Xsq_powers(:,:,2) = X;
normX = norm(X,1);

coef_mat_index = 0;
[B0, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u, coef_mat_index);
normB0 = normest1(double(B0));
normY = normest1(double(Y));

if normX<=s/e
    while normY>normB0*normX^s && s<m 
        B0 = B0 + Y/sfactorial;
        s = s+1;
        normB0 = normest1(double(B0));
        sfactorial = sfactorial*s;
        Xsq_powers(:,:,s) = Y*mp(X);
        Y = Xsq_powers(:,:,s);
        normY = normest1(double(Y));
    end
    r = floor(m/s);
    Bsq_coefmat = zeros(n,n,r+1,'mp');
    Bsq_coefmat(:,:,1) = B0;
    Bsq_coefmat_norm = zeros(1,r+1,'mp');
    usq_prec = zeros(1,r+1,'mp');
    Bsq_coefmat_norm(1) = normB0;
    usq_prec(1) = u;
    for j=1:r
        Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
        Bsq_coefmat_norm(j+1) = normest1(Bsq_coefmat(:,:,j+1)); 
        usq_prec(j+1) = Bsq_coefmat_norm(1)*usq_prec(1)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j));
    end
else
    r = floor(m/s);
    nu = r+1;
    Bsq_coefmat = zeros(n,n,r+1,'mp');
    Bsq_coefmat(:,:,1) = B0;
    Bsq_coefmat_norm = zeros(1,r+1,'mp');
    usq_prec = zeros(1,r+1,'mp');
    Bsq_coefmat_norm(1) = normB0;
    usq_prec(1) = u;
    for j=1:r
        Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
        Bsq_coefmat_norm(j+1) = normest1(Bsq_coefmat(:,:,j+1)); 
        usq_prec(j+1) = max(u, Bsq_coefmat_norm(1)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j))*u);
        if usq_prec(j+1)>=delta*u
            nu = j;
            break;
        end
    end
    for j=nu+1:r
        Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, d_u, j);
        Bsq_coefmat_norm(j+1) = normest1(Bsq_coefmat(:,:,j+1)); 
        usq_prec(j+1) = Bsq_coefmat_norm(nu)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j-nu+1))*u;
    end
end

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

di_vec = double(round(log10(1./usq_prec(2:(r+1))), 0));
Cp = (r*log10(u) + sum(di_vec)) / ((s+r-1)*log10(u)); % the complexity reduction


%% -----------------------       SUBFUNCTIONS       -----------------------
function [q, Xsq_powers, Y] = eval_B0(Xsq_powers, s, d_u, coef_mat_index)
% evaluate the coefficient matrix B_0(X) in precision u

mp.Digits(d_u);
coef_expTaylor = 1./factorial(mp((0:(s-1))+s*coef_mat_index))';

q = coef_expTaylor(1)*Xsq_powers(:,:,1);
for i=1:s-1
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
    q = q + coef_expTaylor(i+1)*Xsq_powers(:,:,i+1);
end
Y = Xsq_powers(:,:,s)*Xsq_powers(:,:,2);

end

function q = eval_coef_mat(Xsq_powers, s, d_u, coef_mat_index)
% form B_i using the powers from Xsq_powers

mp.Digits(d_u);

if coef_mat_index<r
    % form B_i(X), i=1:r-1
    num_term = s-1;
    coef_expTaylor = 1./factorial(mp((0:num_term)+s*coef_mat_index))';
else
    % form B_r(X)
    num_term = m-s*r;
    coef_expTaylor = 1./factorial(mp((0:num_term)+s*coef_mat_index))';
end
reshaped_vector = reshape(coef_expTaylor,[1,1,num_term+1]);
q = sum(bsxfun(@times, Xsq_powers(:,:,1:(num_term+1)),reshaped_vector),3);
end
end