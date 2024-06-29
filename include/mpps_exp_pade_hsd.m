function [p, s, cr, usq_prec] = mpps_exp_pade_hsd(X, m, padepoly)
%PS_EXP_PADE_HSD    The mixed-precision paterson-stockmeyer algorithm 
% for the Pade approximants of the matrix exponential, it uses only 
% possibly double, single, and half precisions.
%
% The degree of the polynomials in the numerator and denominator of Pade 
% approximant is specified by m, PADEPOLY specifies the kind by
% PADEPOLY = NUMERATOR or PADEPOLY = DENOMINATOR.
% P is the computed polynomial, the parameter in the PS scheme is 
% S = ceil(sqrt(m)), USQ_PREC is a vector that stores the used
% precisions, and CR is the complexity reduction.

n = size(X,1); % size of the matrix
s = ceil(sqrt(m));  % the (default) parameter s in the PS scheme
delta = 10; % precision switching parameter

% vector of Pade coefficients
c_pade_coefficients = get_pade_coefficients(m, padepoly);

u = eps('double')/2;
u_sgl = eps('single')/2;
u_hlf = 2^(-8); % bfloat16
options.format = 'b'; options.round = 1; options.subnormal = 1; % bfloat16
chop([],options)

% u_hlf = 2^(-11); % issue of excessive underflow for fp16
% options.format = 'h'; options.round = 1; options.subnormal = 1; % fp16 


Xsq_powers = zeros(n,n,s);
Xsq_powers(:,:,1) = eye(n);
Xsq_powers(:,:,2) = X;

[B0, Xsq_powers, Y] = eval_B0(Xsq_powers, s);
normB0 = normest1(B0);
normY = normest1(Y);

r = floor(m/s);
nu = r+1;
Bsq_coefmat = zeros(n,n,r+1);
Bsq_coefmat(:,:,1) = B0;
Bsq_coefmat_norm = zeros(1,r+1);
usq_prec = zeros(1,r+1);
Bsq_coefmat_norm(1) = normB0;
usq_prec(1) = u;
for j=1:r
    Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, j);
    Bsq_coefmat_norm(j+1) = normest1(single(Bsq_coefmat(:,:,j+1))); 
    usq_prec(j+1) = max(u, Bsq_coefmat_norm(1)/(Bsq_coefmat_norm(j+1)*mp(normY)^(j))*u);
    if usq_prec(j+1)>=eps('single')/2/delta
        nu = j; 
        break;
    end
end
for j=nu+1:r
        Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, j);
        Bsq_coefmat_norm(j+1) = normest1(Bsq_coefmat(:,:,j+1)); 
        usq_prec(j+1) = Bsq_coefmat_norm(nu)/(Bsq_coefmat_norm(j+1)*normY^(j-nu+1))*u;
end

% compute the polynomial p using Horner's method
p = Bsq_coefmat(:,:,r+1);
for k=r:-1:1
    p = matmul(p, Y, usq_prec(k+1)); % matrix multiply in double/single/half precision
    p = matadd(p, Bsq_coefmat(:,:,k), usq_prec(k)); % matrix addition in double/single/half precision
end

cr = comput_complex_reduc_hsd(usq_prec);

%% -----------------------       SUBFUNCTIONS       -----------------------
function [q, Xsq_powers, Y] = eval_B0(Xsq_powers, s)
% evaluate the coefficient matrix B_0(X) in double precision


coef_expPade = c_pade_coefficients((0:(s-1))+1)';

q = coef_expPade(1)*Xsq_powers(:,:,1);
for i=1:s-1
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
    q = q + coef_expPade(i+1)*Xsq_powers(:,:,i+1);
end
Y = Xsq_powers(:,:,s)*Xsq_powers(:,:,2);

end

function [q, Xsq_powers] = eval_coef_mat(Xsq_powers, s, coef_mat_index)
% evaluate B_i using Xsq_powers, no O(n^3) operations so this is done 
% completely in double precision

if coef_mat_index<r
    % form B_i(X), i=1:r-1
    num_term = s-1;
    coef_expPade =  c_pade_coefficients((0:num_term)+s*coef_mat_index+1)';
else
    % form B_r(X)
    num_term = m-s*r;
    coef_expPade = c_pade_coefficients((0:num_term)+s*coef_mat_index+1)';
end
reshaped_vector = reshape(coef_expPade,[1,1,num_term+1]);
q = sum(bsxfun(@times, Xsq_powers(:,:,1:(num_term+1)),reshaped_vector),3);

end

function C = matmul(A, B, ui)
% matrix multiplication in double/single/half precision

if ui<u_sgl % use double precision
        C = A*B;
elseif ui<u_hlf % use single precision
        C = double(single(A)*single(B));
else % use half precision
        C = matmul_half(A, B);
end
end

function C = matadd(A, B, ui)
% matrix addition in double/single/half precision

if ui<u_sgl % use double precision
    C = A+B;
elseif ui<u_hlf % use single precision
    C = double(single(A)+single(B));
else % use half precision
    C = matadd_half(A, B);
end
end

function C = matmul_half(A, B)
% matrix multiplication in half precision
if isreal(A) && isreal(B) % both A and B are real matrices
    C = matmul_half_real(A, B);
else
    imagA = imag(A);
    realA = real(A); 
    imagB = imag(B);
    realB = real(B); 
    C = matmul_half_real(realA, realB) - matmul_half_real(imagA, imagB) + ...
        1i*(matmul_half_real(realA, imagB) + matmul_half_real(imagA, realB));
end
end

function C = matadd_half(A, B)
% matrix addition in half precision
if isreal(A) && isreal(B) % both A and B are real matrices
    C = matadd_half_real(A, B);
else
    imagA = imag(A);
    realA = real(A); 
    imagB = imag(B);
    realB = real(B); 
    C = matadd_half_real(realA, realB) + 1i*matadd_half_real(imagA, imagB);
end
end 

function C = matmul_half_real(A, B)
% matrix multiplication between real matrices in half precision
    A = chop(A);
    B = chop(B);
    C = zeros(n,'double');
    for i = 1:n
        C = chop(C + chop(A(:,i)*B(i,:)));
    end
end

function C = matadd_half_real(A, B)
% matrix addition between real matrices in half precision
    C = chop(A+chop(B));
end

function c = get_pade_coefficients(m, options)
%get_pade_coefficients Coefficients of Pade approximant
%    C = get_pade_coefficients returns coefficients of
%    numerator/denominator of [m/m] Pade approximant, where m = 3,5,7,9,13.

if isequal(options,'numerator')
    switch m
    case 3
        c = [1, 0.5, 0.1, 0.008333333333333333];
    case 5
        c = [1, 0.5, 0.1111111111111111, 0.01388888888888889, 0.000992063492063492, 3.306878306878307e-05];
    case 7
        c = [1, 0.5, 0.1153846153846154, 0.01602564102564102, 0.001456876456876457, 8.741258741258741e-05...
            3.237503237503237e-06, 5.781255781255781e-08];
    case 9
        c = [1, 0.5, 0.1176470588235294, 0.01715686274509804, 0.001715686274509804, 0.0001225490196078431...
            6.28456510809452e-06, 2.244487538605186e-07, 5.101108042284513e-09, 5.66789782476057e-11];
    case 13
        c = [1, 0.5, 0.12, 0.01833333333333333, 0.001992753623188406, 0.0001630434782608696,...
            1.0351966873706e-05, 5.175983436853002e-07, 2.043151356652501e-08, 6.306022705717595e-10...
            1.48377004840414e-11, 2.529153491597966e-13, 2.810170546219962e-15, 1.544049750670309e-17];
    end
end

if isequal(options,'denominator')
    switch m
    case 3
        c = [1, -0.5, 0.1, -0.008333333333333333];
    case 5
        c = [1, -0.5, 0.1111111111111111, -0.01388888888888889, 0.000992063492063492, -3.306878306878307e-05];
    case 7
        c = [1, -0.5, 0.1153846153846154, -0.01602564102564102, 0.001456876456876457, -8.741258741258741e-05...
            3.237503237503237e-06, -5.781255781255781e-08];
    case 9
        c = [1, -0.5, 0.1176470588235294, -0.01715686274509804, 0.001715686274509804, -0.0001225490196078431...
            6.28456510809452e-06, -2.244487538605186e-07, 5.101108042284513e-09, -5.66789782476057e-11];
    case 13
        c = [1, -0.5, 0.12, -0.01833333333333333, 0.001992753623188406, -0.0001630434782608696, 1.0351966873706e-05...
            -5.175983436853002e-07, 2.043151356652501e-08, -6.306022705717595e-10, 1.48377004840414e-11, -2.529153491597966e-13...
            2.810170546219962e-15, -1.544049750670309e-17];
    end
end

end

function cr = comput_complex_reduc_hsd(usq_prec)
% compute complexity reduction in percentage from the vector USQ_PREC which
% stores the used precisions.

num_half = sum(usq_prec>=u_hlf);
num_sgl = sum(usq_prec>=u_sgl) - num_half;
num_dbl = length(usq_prec) - 1 - num_sgl - num_half;

cr = 1 - (4*(s-1+num_dbl)+2*num_sgl+num_half) / (4*(s-1+r));

end
end