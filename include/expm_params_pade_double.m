function [s, m] = expm_params_pade_double(A)
%EXPM_PARAMS_PADE_DOUBLE Obtain scaling parameter and order of the Pade 
% approximant.
%
%   taken from:
%   A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm
%   for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3),
%   (2009), pp. 970-989.

% Coefficients of backward error function.
coeff = [1/100800, 1/10059033600, 1/4487938430976000,...
     1/5914384781877411840000, 1/113250775606021113483283660800000000];

s = 0;
% m_val is one of [3 5 7 9 13];
% theta_m for m=1:13.
theta = [%3.650024139523051e-008
    %5.317232856892575e-004
    1.495585217958292e-002  % m_vals = 3
    %8.536352760102745e-002
    2.539398330063230e-001  % m_vals = 5
    %5.414660951208968e-001
    9.504178996162932e-001  % m_vals = 7
    %1.473163964234804e+000
    2.097847961257068e+000  % m_vals = 9
    %2.811644121620263e+000
    %3.602330066265032e+000
    %4.458935413036850e+000
    5.371920351148152e+000];% m_vals = 13

Apowers{2} = A*A;
Apowers{4} = Apowers{2}*Apowers{2};
Apowers{6} = Apowers{2}*Apowers{4};
d4 = norm(Apowers{4},1)^(1/4);
d6 = norm(Apowers{6},1)^(1/6);
eta1 = max(d4, d6);
if (eta1 <= theta(1) && ell(A, coeff(1), 3) == 0)
    m = 3;
    return;
end
if (eta1 <= theta(2) && ell(A, coeff(2), 5) == 0)
    m = 5;
    return;
end

isSmall = size(A,1) < 150; %Compute matrix power explicitly
if isSmall
    d8 = norm(Apowers{4}*Apowers{4},1)^(1/8);
else
    d8 = normAm(Apowers{4}, 2)^(1/8);
end
eta3 = max(d6, d8);
if (eta3 <= theta(3) && ell(A, coeff(3), 7) == 0)
    m = 7;
    return;
end
if (eta3 <= theta(4) && ell(A, coeff(4), 9) == 0)
    m = 9;
    return;
end
if isSmall
    d10 = norm(Apowers{4}*Apowers{6},1)^(1/10);
else
    d10 = normAm(Apowers{2}, 5)^(1/10);
end
eta4 = max(d8, d10);
eta5 = min(eta3, eta4);
s = max(ceil(log2(eta5/theta(5))), 0);
s = s + ell(A/2^s, coeff(5), 13);
if isinf(s)
    % Overflow in ell subfunction. Revert to old estimate.
    [t, s] = log2(norm(A,1)/theta(end));
    s = s - (t == 0.5); % adjust s if normA/theta(end) is a power of 2.
end
m = 13;
end

%% Subfunctions
function t = ell(A, coeff, m_val)
%ell Function needed to compute optimal parameters.
scaledT = coeff.^(1/(2*m_val+1)) .* abs(A);
alpha = normAm(scaledT,2*m_val+1)/norm(A,1);
t = max(ceil(log2(2*alpha/eps(class(alpha)))/(2*m_val)),0);
end

function [c,mv] = normAm(A,m)
%NORMAM   Estimate of 1-norm of power of matrix.
%   NORMAM(A,m) estimates norm(A^m,1).
%   If A has nonnegative elements the estimate is exact.
%   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
%   matrix-vector products computed involving A or A^*.

%   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
%   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
%   970-989, 2009.

%   Awad H. Al-Mohy and Nicholas J. Higham, April 19, 2010.

n = length(A);
if isequal(A,abs(A))
    e = ones(n,1);
    for j=1:m         % for positive matrices only
        e = A'*e;
    end
    c = norm(e,inf);
    mv = m;
else
    [c,v,w,it] = normest1(@afun_power);
    mv = it(2)*2*m; % Since t = 2.
end

  function Z = afun_power(flag,X)
       %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

       if isequal(flag,'dim')
          Z = n;
       elseif isequal(flag,'real')
          Z = isreal(A);
       else

          [p,q] = size(X);
          if p ~= n, error('Dimension mismatch'), end

          if isequal(flag,'notransp')
             for i = 1:m, X = A*X; end
          elseif isequal(flag,'transp')
             for i = 1:m, X = A'*X; end
          end

          Z = X;
       end
  end
end