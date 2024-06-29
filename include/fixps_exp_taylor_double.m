function p = fixps_exp_taylor_double(X, m)
%FIXPS_EXP_TAYLOR_DOUBLE    The fixed-precision Paterson-Stockmeyer scheme 
% for the Taylor approximants of the matrix exponential in double precision.
% 
% The degree of the Taylor approximant is specified by M, and the computed 
% polynomial is returned as P.


n = size(X,1);
s = ceil(sqrt(m));
r = floor(m/s);

Xsq_powers = zeros(n,n,s,'double');
Xsq_powers(1:n+1:n^2) = 1;
Xsq_powers(:,:,2) = X;
I = Xsq_powers(:,:,1);

c = 1./factorial((0:m))';


% Compute first s+1 powers.
for i=2:s
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
end
     
% Evaluate last polynomial of degree m mod s.
B = c(m+1) * Xsq_powers(:,:,m-s*r+1);
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
    B = zeros(n, 'double');
    B = B + c(s*kk+1) * I;
    for j=1:s-1
      B = B + c(s*kk+j+1) * Xsq_powers(:,:,j+1);
    end
    p = p * Xsq_powers(:,:,s+1) + B;
end
end