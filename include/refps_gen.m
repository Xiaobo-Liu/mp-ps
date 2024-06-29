function p = refps_gen(X, m, coef, digs_wkprec)
%REFPS_GEN   Compute the reference solution of a general polynomial of  
% the matrix X.
%
% The input COEF is a (m+1)-vector which stores the coefficients of the
% matrix polynomial up to degree M and, if the computing environment is 
% arbitrary precision, the equivalent decimal digits of the working 
% precision should be pass as DIGS_WKPREC.
% The computed polynomial is returned as P.


digs_ref = 2 * digs_wkprec;
mp.Digits(digs_ref);

s = ceil(sqrt(m));
r = floor(m/s);

n = size(X,1);
X = mp(X);

Xsq_powers = zeros(n,n,s,'mp');
Xsq_powers(1:n+1:n^2) = 1;
Xsq_powers(:,:,2) = X;
I = Xsq_powers(:,:,1);


% Compute first s+1 powers.
for i=2:s
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
end
     
% Evaluate last polynomial of degree m mod s.
B = mp(coef(m+1) * Xsq_powers(:,:,m-s*r+1));
for j=m-1:-1:s*r
   if j == s*r
       B = B + coef(s*r+1)*I;
   else      
       B = B + coef(j+1) * Xsq_powers(:,:,m-s*r-(m-j)+1);
   end
end

% Evaluate polynomials of degree s-1 and evaluate main polynomial using
% Horner's method.
p = B;
for kk=r-1:-1:0
    B = zeros(n, 'mp');
    B = B + coef(s*kk+1) * I;
    for j=1:s-1
      B = B + coef(s*kk+j+1) * Xsq_powers(:,:,j+1);
    end
    p = p * Xsq_powers(:,:,s+1) + B;
end
end