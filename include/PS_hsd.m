function [p, s, usq_prec] = PS_hsd(X, m, delta)
% The mixed-precision paterson-stockmeyer algorithm which uses only double,
% single, and half precisions.

% The degree of the Taylor approximant is specified by m, 
% the working precision is double,
% and delta = 10 is recommended.

n = size(X,1); % size of the matrix
s = ceil(sqrt(m));  % the (default) parameter s in the PS scheme
sfactorial = factorial(s);
e = exp(1);

u = eps('double')/2;
u_sgl = eps('single')/2;

% u_hlf = 2^(-11); % fp16
% options.format = 'h'; options.round = 1; options.subnormal = 1; % fp16

u_hlf = 2^(-8); % bfloat16
options.format = 'b'; options.round = 1; options.subnormal = 1; % bfloat16

chop([],options)

Xsq_powers = zeros(n,n,s);
Xsq_powers(:,:,1) = eye(n);
Xsq_powers(:,:,2) = X;
normX = norm(X,1);

coef_mat_index = 0;
[B0, Xsq_powers, Y] = eval_B0(Xsq_powers, s, coef_mat_index);
normB0 = normest1(B0);
normY = normest1(Y);

if normX<=s/e
    while normY>normB0*normX^s && s<m 
        B0 = B0 + Y/sfactorial;
        s = s+1;
        normB0 = normest1(B0);
        sfactorial = sfactorial*s;
        Xsq_powers(:,:,s) = Y*X;
        Y = Xsq_powers(:,:,s);
        normY = normest1(Y);
    end
    r = floor(m/s);
    Bsq_coefmat = zeros(n,n,r+1);
    Bsq_coefmat(:,:,1) = B0;
    Bsq_coefmat_norm = zeros(1,r+1);
    usq_prec = zeros(1,r+1);
    Bsq_coefmat_norm(1) = normB0;
    usq_prec(1) = u;
    for j=r:-1:1
        Bsq_coefmat(:,:,j+1) = eval_coef_mat(Xsq_powers, s, j);
        Bsq_coefmat_norm(j+1) = normest1(single(Bsq_coefmat(:,:,j+1))); 
        usq_prec(j+1) = Bsq_coefmat_norm(1)/(Bsq_coefmat_norm(j+1)*normY^(j))*usq_prec(1);
    end
else
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
end

% compute the polynomial p using Horner's method
p = Bsq_coefmat(:,:,r+1);
for k=r:-1:1
    p = matmul(p, Y, usq_prec(k+1)); % matrix multiply in double/single/half precision
    p = matadd(p, Bsq_coefmat(:,:,k), usq_prec(k)); % matrix addition in double/single/half precision
end


%% -----------------------       SUBFUNCTIONS       -----------------------
function [q, Xsq_powers, Y] = eval_B0(Xsq_powers, s, coef_mat_index)
% evaluate the coefficient matrix B_0(X) in double precision


coef_expTaylor = 1./factorial((0:(s-1)+s*coef_mat_index))';

q = coef_expTaylor(1)*Xsq_powers(:,:,1);
for i=1:s-1
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
    q = q + coef_expTaylor(i+1)*Xsq_powers(:,:,i+1);
end
Y = Xsq_powers(:,:,s)*Xsq_powers(:,:,2);

end

function [q, Xsq_powers] = eval_coef_mat(Xsq_powers, s, coef_mat_index)
% evaluate B_i using Xsq_powers, no O(n^3) operations so this is done 
% completely in double precision

if coef_mat_index<r
    % form B_i(X), i=1:r-1
    num_term = s-1;
    coef_expTaylor = 1./factorial((0:num_term)+s*coef_mat_index)';
else
    % form B_r(X)
    num_term = m-s*r;
    coef_expTaylor = 1./factorial((0:num_term)+s*coef_mat_index)';
end
reshaped_vector = reshape(coef_expTaylor,[1,1,num_term+1]);
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

end