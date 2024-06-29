function p = refps_exp_taylor(X, m, digs_wkprec)
%REFPS_EXP_TAYLOR   Compute the reference solution of a Taylor approximant 
% of the matrix exponential at the matrix X.
%
% The degree of the Taylor approximant is specified by M and, if the 
% computing environment is arbitrary precision, the equivalent decimal
% digits of the working precision should be pass as DIGS_WKPREC.
% The computed polynomial is returned as P.

if nargin < 3 || isempty(digs_wkprec) % the working precision is double
    digs_ref = 100; 
else 
    digs_ref = 2 * digs_wkprec;
end

mp.Digits(digs_ref);
n = size(X,1);
X = mp(X);

s = ceil(sqrt(m));
r = floor(m/s);

Xsq_powers = zeros(n,n,s,'mp');
Xsq_powers(1:n+1:n^2) = 1;
Xsq_powers(:,:,2) = X;
I = Xsq_powers(:,:,1);

c = 1./factorial(mp((0:m)))';

% Compute first s+1 powers.
for i=2:s
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
end
     
% Evaluate last polynomial of degree m mod s.
B = mp(c(m+1) * Xsq_powers(:,:,m-s*r+1));
for j=m-1:-1:s*r
   if j == s*r
       B = B + c(s*r+1)*I;
   else      
       B = B + c(j+1) * Xsq_powers(:,:,m-s*r-(m-j)+1);
   end
end

% Evaluate polynomials of degree s-1 and evaluate main polynomial using
% Horner's method.
p = B;
for kk=r-1:-1:0
    B = zeros(n, 'mp');
    B = B + c(s*kk+1) * I;
    for j=1:s-1
      B = B + c(s*kk+j+1) * Xsq_powers(:,:,j+1);
    end
    p = p * Xsq_powers(:,:,s+1) + B;
end
end