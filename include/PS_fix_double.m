function F = PS_fix_double(X, mm)
% The fixed-precision PS scheme for the Taylor approximant of degree m of 
% exp at the matrix X in double precision.



n = size(X,1);
ss = ceil(sqrt(mm));
rr = floor(mm/ss);

Xsq_powers = zeros(n,n,ss,'double');
Xsq_powers(:,:,1) = eye(n,'double');
Xsq_powers(:,:,2) = X;
I = Xsq_powers(:,:,1);

c = 1./factorial((0:mm))';


% Compute first s+1 powers.
for i=2:ss
    Xsq_powers(:,:,i+1) = Xsq_powers(:,:,i)*Xsq_powers(:,:,2);
end
     
% Evaluate last polynomial of degree m mod s.
B = c(mm+1) * Xsq_powers(:,:,mm-ss*rr+1);
for j=mm-1:-1:ss*rr
   if j == ss*rr
       B = B + c(ss*rr+1)*I;
   else      
       B = B + c(j+1) * Xsq_powers(:,:,mm-ss*rr-(mm-j)+1);
   end
end

% Evaluate polynomials of degree s-1 and evaluate main polynomial using
% Horner's method.
F = B;
for kk=rr-1:-1:0
    B = zeros(n, 'double');
    B = B + c(ss*kk+1) * I;
    for j=1:ss-1
      B = B + c(ss*kk+j+1) * Xsq_powers(:,:,j+1);
    end
    F = F * Xsq_powers(:,:,ss+1) + B;
end
end