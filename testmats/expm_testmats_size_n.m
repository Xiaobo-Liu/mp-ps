function [A, n_mats] = expm_testmats_size_n(k, n)
%EXPM_TESTMATS_VARN  Test matrices for matrix exponential with variable 
%size (size(A,1) == n). This is a special subset of our EXPM_TESTMATS.
%
%   [A, NMATS] = EXPM_TESTMATS_VARN(K,N) selects the K'th test matrix.
%   NMATS is the number of test matrices available.
%   N sets the dimension of the matrices if the size is undetermined.
%
%   In total there are 23 + 17 + 29 + 4 = 73 test matrices.


    indices = [101:123, 201:217, 301:329, 401:404];

    n_mats = length(indices);
    if nargin < 1
        A = [];
        return;
    end
    if nargin < 2
        error('The dimension n of test matrices is not specified.');
    end
    local_index = indices(k);


    switch local_index

        %% Positive real eigenvalues
      case 101, sA='cauchy';A=gallery(sA,n); %symm % optional parameters
      case 102, sA='condex';A=gallery(sA,n,3); %lower triangular
      case 103, sA='condex';A=gallery(sA,n,4,100); %symm real
      case 104, sA='dorr';A=full(gallery(sA,n,100));
      case 105, sA='dramadah';A=gallery(sA,n,2);
      case 106, sA='frank';A=gallery(sA,n);
      case 107, sA='gcdmat';A=gallery(sA,n); %symm real
      case 108, sA='grcar';A=gallery(sA,n);
      case 109, sA='hanowa';A=gallery(sA,n);
      case 110, A=hilb(n); % sA='hilb'; symm real

      case 111, sA='invhess';A=gallery(sA,n);
      case 112, sA='jordbloc';A=gallery(sA,n,1); %symm real
      case 113, sA='kahan';A=gallery(sA,n);
      case 114, sA='lehmer';A=gallery(sA,n); %symm real
      case 115, sA='minij';A=gallery(sA,n); %symm real
      case 116, sA='moler';A=gallery(sA,n); %symm real
      case 117, sA='parter';A=gallery(sA,n);
      case 118, sA='pei';A=gallery(sA,n);
      case 119, sA='prolate';A=gallery(sA,n,1); %symm real Toeplitz
      case 120, sA='randcorr';A=gallery(sA,n); %symm real

      case 121, sA='sampling';A=gallery(sA,n);
      case 122, sA='toeppd';A=gallery(sA,n); %symm real
      case 123, sA='tridiag';A=full(gallery(sA,n));

        %% Real eigenvalues
      case 201, sA='fiedler';A=gallery(sA,n);
      case 202, sA = 'house'; [v, beta] = gallery(sA, n);
        A = eye(n) - beta * (v * v');
      case 203, sA='jordbloc';A=gallery(sA,n,2); %symm real
      case 204, sA='kms';A=gallery(sA,n); %symm real
      case 205, sA='lesp';A=gallery(sA,n);
      case 206, sA='lotkin';A=gallery(sA,n);
      case 207, sA='orthog';A=gallery(sA,n,1); %symm real
      case 208, sA='orthog';A=gallery(sA,n,2); %symm real
      case 209, sA='orthog';A=gallery(sA,n,5); %symm real
      case 210, sA='orthog';A=gallery(sA,n,6); %symm real

      case 211, sA='orthog';A=gallery(sA,n,-1); %symm real
      case 212, sA='redheff';A=gallery(sA,n);
      case 213, sA='riemann';A=gallery(sA,n);
      case 214, sA='ris';A=gallery(sA,n,1e1);
      case 215, sA='clement';A=gallery(sA,n,1); %symm
      case 216, A = gallery('triw',n,1,1); 
          A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2); % \cite[p.~10, Ex III]{pang85}
      case 217, alpha = 1; beta = 1; % \cite{kuda10}, no values are given in the paper.
          A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
          A(1,2) = beta; A(n,n-1) = beta;  % very ill conditioned
      % case 218, sA='binomial'; A=gallery(sA,n); % causing overflow with n=500, removed

        %% Complex eigenvalues
      case 301, sA='chebspec';A=gallery(sA,n);
      case 302, sA='chebvand';A=gallery(sA,n);
      case 303, sA='chow';A=gallery(sA,n);
      case 304, sA='circul';A=gallery(sA,n);
      case 305, sA='cycol';A=gallery(sA,n);
      case 306, sA='dramadah';A=gallery(sA,n,1);
      case 307, sA='dramadah';A=gallery(sA,n,3);
      case 308, sA='forsythe';A=gallery(sA,n);
      case 309, sA='leslie'; A = gallery(sA,n);
      case 310, sA='leslie';A=gallery(sA,n);

      case 311, sA='normaldata';A=gallery(sA,n,10);
      case 312, sA='orthog';A=gallery(sA,n,3); %symm non herm
      case 313, sA='orthog';A=gallery(sA,n,4);
      case 314, sA='orthog';A=gallery(sA,n,-2);
      case 315, sA='randcolu';A=gallery(sA,n);
      case 316, sA='randhess';A=gallery(sA,n);
      case 317, sA='rando';A=gallery(sA,n,1);
      case 318, sA='rando';A=gallery(sA,n,2);
      case 319, sA='rando';A=gallery(sA,n,3);
      case 320, sA='randsvd';A=gallery(sA,n,1);

      case 321, sA='randsvd';A=gallery(sA,n,2);
      case 322, sA='randsvd';A=gallery(sA,n,3);
      case 323, sA='randsvd';A=gallery(sA,n,4);
      case 324, sA='randsvd';A=gallery(sA,n,5);
      case 325, sA='smoke';A=gallery(sA,n);
      case 326, sA='smoke';A=gallery(sA,n,1);
      case 327, sA='toeppen';A=full(gallery(sA,n));
      case 328, sA='uniformdata';A=gallery(sA,n,1000);
      case 329, A = gallery('triw',n,1);  m = (n-1)/2; % \cite[p.~9, Ex II]{pang85}
          A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1); % Fasi's interpretation of their matrix for arbitrary n.
          for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end
      % case 330, sA='invol';A=gallery(sA,n); % causing overflow with n=500, removed

        %% Singular matrices
      case 401, sA='gearmat';A=gallery(sA,n);
      case 402, A = 1e2*triu(randn(n),1); % stricly uppertriangular so nilpotent.
      case 403, A = zeros(n); A(n+1:n+1:n^2) = 1:n-1; % log of Cholesky factor of Pascal matrix. \cite{edst03}.
      case 404, v = eye(n,1); A = eye(n) - v*v'; % identity matrix with (1,1) entry replaced by 0
    end   
end