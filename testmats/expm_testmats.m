function [A, n_mats] = expm_testmats(k, n)
%EXPM_TESTMATS Test matrices for the matrix exponential.
%
%   [A, NMATS] = EXPM_TESTMATS(K,N) selects the K'th test matrix.
%   NMATS is the number of test matrices available. N sets the dimension of 
%   the matrices if the size is undetermined; by default N=20.
%   In total, there are 79+42 = 121 test matrices, including all matrices
%   from EXPM_TESTMATS_VARDIM and all matrices from the MATRICES-EXPM 
%   collection at 
%   https://github.com/Xiaobo-Liu/matrices-expm
%
%   Matrices from the MATRICES-EXPM collection are in Group 5, in addition
%   to the other four groups inherited from EXPM_TESTMATS_VARDIM. 
%   Groups 1-4 consist of matrices of variable size that are more general 
%   and are from the Anymatrix MATLAB toolbox
%   https://github.com/mmikaitis/anymatrix
%   and they are classified by having
%   1) positive real eigenvalues, 2) real eigenvalues, 
%   3) complex eigenvalues, and 4) zero as an eigenvalue.
%   Matrices from the MATRICES-EXPM collection are in group 5.

    indices = [101:126, 201:217, 301:333, 401:403, 501:542];

    n_mats = length(indices);
  
    if nargin < 1
        A = [];
        return;
    end
    if nargin < 2
        n = 20; % the default dimention if it is not specified
    end
    local_index = indices(k);

    switch local_index

        % Positive real eigenvalues
      case 101, sA='core/stoch_cesaro'; A=anymatrix(sA,n); 
      case 102, sA='core/triminsval01'; A=anymatrix(sA,n);
      case 103, sA='gallery/cauchy'; A=anymatrix(sA,n); %symm
      case 104, sA='gallery/condex'; A=anymatrix(sA,n,3); %mode 3, lower triangular
      case 105, sA='gallery/condex'; A=anymatrix(sA,n,4,100); %mode 4, symm real
      case 106, sA='gallery/dorr';   A=full(anymatrix(sA,n,100)); %theta=100 to make elements big
      case 107, sA='gallery/dramadah'; A=anymatrix(sA,n,2); %mode 2, upper triangular and Toeplitz
      case 108, sA='gallery/frank'; A=anymatrix(sA,n);
      case 109, sA='gallery/gcdmat'; A=anymatrix(sA,n); %symm real
      case 110, sA='gallery/grcar'; A=anymatrix(sA,n);
      case 111, sA='gallery/hanowa'; A=anymatrix(sA,n);
      case 112, sA='gallery/invhess'; A=anymatrix(sA,n);
      case 113, sA='gallery/jordbloc';A=anymatrix(sA,n); 
      case 114, sA='gallery/kahan'; A=anymatrix(sA,n);
      case 115, sA='gallery/lehmer'; A=anymatrix(sA,n); %symm real
      case 116, sA='gallery/minij'; A=anymatrix(sA,n); %symm real
      case 117, sA='gallery/moler'; A=anymatrix(sA,n); %symm real
      case 118, sA='gallery/parter'; A=anymatrix(sA,n);
      case 119, sA='gallery/pei'; A=anymatrix(sA,n);
      case 120, sA='gallery/prolate'; A=anymatrix(sA,n); %symm real Toeplitz
      case 121, sA='gallery/randcorr'; A=anymatrix(sA,n); %symm real
      case 122, sA='gallery/sampling'; A=anymatrix(sA,n); 
      case 123, sA='gallery/toeppd'; A=anymatrix(sA,n); %symm real
      case 124, sA='gallery/tridiag'; A=full(anymatrix(sA,n)); %symm real M-matrix
      case 125, sA='gallery/triw'; A=anymatrix(sA,n); 
      case 126, sA='matlab/hilb'; A=anymatrix(sA,n); %symm real

        % Real eigenvalues
      case 201, sA='core/blockhouse'; A=anymatrix(sA,n); %block Householder matrix
      case 202, sA='core/stoch_revtri'; A=anymatrix(sA,n);
      case 203, sA='gallery/clement'; A=anymatrix(sA,n); %nonsymm
      case 204, sA='gallery/clement'; A=anymatrix(sA,n,1); %k=1, symm
      case 205, sA='gallery/fiedler'; A=anymatrix(sA,n); %symm real
      case 206, sA='gallery/kms'; A=anymatrix(sA,n); %symm real
      case 207, sA='gallery/lesp'; A=anymatrix(sA,n);
      case 208, sA='gallery/lotkin'; A=anymatrix(sA,n);
      case 209, sA='gallery/orthog'; A=anymatrix(sA,n,1);%type 1, symm real
      case 210, sA='gallery/orthog'; A=anymatrix(sA,n,2); %symm real
      case 211, sA='gallery/orthog'; A=anymatrix(sA,n,3); %symm real
      case 212, sA='gallery/orthog'; A=anymatrix(sA,n,5); %symm real
      case 213, sA='gallery/orthog'; A=anymatrix(sA,n,6); %symm real
      case 214, sA='gallery/orthog'; A=anymatrix(sA,n,-1); %symm real
      case 215, sA='gallery/redheff'; A=anymatrix(sA,n);
      case 216, sA='gallery/riemann'; A=anymatrix(sA,n);
      case 217, sA='gallery/ris'; A=anymatrix(sA,n); %symm real
     
        % Complex eigenvalues
      case 301, sA='core/hess_orth'; A=anymatrix(sA,n);
      case 302, sA='core/orthog_cauchy'; A=anymatrix(sA,n); %Orthogonal Cauchy-like
      case 303, sA='core/pick'; A=anymatrix(sA,n);
      case 304, sA='core/soules'; A=anymatrix(sA,n);
      case 305, sA='core/stoch_compan'; A=anymatrix(sA,n);
      case 306, sA='core/tournament'; A=anymatrix(sA,n);
      case 307, sA='gallery/chebspec'; A=anymatrix(sA,n);
      case 308, sA='gallery/chebvand'; A=anymatrix(sA,n);
      case 309, sA='gallery/chow'; A=anymatrix(sA,n);
      case 310, sA='gallery/circul'; A=anymatrix(sA,n);
      case 311, sA='gallery/cycol'; A=anymatrix(sA,n);
      case 312, sA='gallery/dramadah'; A=anymatrix(sA,n);
      case 313, sA='gallery/dramadah'; A=anymatrix(sA,n,3); %interesting e'val distri
      case 314, sA='gallery/forsythe'; A=anymatrix(sA,n);
      case 315, sA='gallery/leslie'; A=anymatrix(sA,n);
      case 316, sA='gallery/normaldata'; A=anymatrix(sA,n,10); %random seed=10
      case 317, sA='gallery/orthog'; A=anymatrix(sA,n,4); %type 4
      case 318, sA='gallery/orthog'; A=anymatrix(sA,n,-2);%type -2  
      case 319, sA='gallery/randcolu'; A=anymatrix(sA,n);
      case 320, sA='gallery/randhess';A=anymatrix(sA,n);
      case 321, sA='gallery/randjorth'; A=anymatrix(sA,ceil(n/2),n-ceil(n/2),1e4); %to make size n
      case 322, sA='gallery/rando'; A=anymatrix(sA,n,1); %mode 1
      case 323, sA='gallery/rando'; A=anymatrix(sA,n,2);
      case 324, sA='gallery/rando'; A=anymatrix(sA,n,3);
      case 325, sA='gallery/randsvd'; A=anymatrix(sA,n,1e6,1); %kappa=1e6; mode 1
      case 326, sA='gallery/randsvd'; A=anymatrix(sA,n,1e6,2); %kappa=1e6; mode 2
      case 327, sA='gallery/randsvd'; A=anymatrix(sA,n,1e6,3); %kappa=1e6; mode 3
      case 328, sA='gallery/randsvd'; A=anymatrix(sA,n,1e6,4); %kappa=1e6; mode 4
      case 329, sA='gallery/randsvd'; A=anymatrix(sA,n,1e6,5); %random kappa
      case 330, sA='gallery/smoke'; A=anymatrix(sA,n);
      case 331, sA='gallery/smoke'; A=anymatrix(sA,n,1); %A(n,1)=zero.
      case 332, sA='gallery/toeppen'; A=full(anymatrix(sA,n)); %pentadiagonal Toeplitz 
      case 333, sA='gallery/uniformdata'; A=anymatrix(sA,n,10); %random seed=10
   
        % Zero eigenvalues
      case 401, sA='core/creation'; A=anymatrix(sA,n); %nilpotent
      case 402, sA='core/hessfull01'; A=anymatrix(sA,n); %rank=n-1
      case 403, sA='gallery/gearmat'; A=anymatrix(sA,n);
      
        % Matrices from matrices-expm collection
      case 501, sA='expm/alhi09r1'; A=anymatrix(sA);
      case 502, sA='expm/alhi09r2'; A=anymatrix(sA);
      case 503, sA='expm/alhi09r3'; A=anymatrix(sA);
      case 504, sA='expm/alhi09r4'; A=anymatrix(sA);
      case 505, sA='expm/dahi03'; A=anymatrix(sA);
      case 506, sA='expm/dipa00'; A=anymatrix(sA);
      case 507, sA='expm/edst04'; A=anymatrix(sA);
      case 508, sA='expm/eigt7'; A=anymatrix(sA);
      case 509, sA='expm/fahi19r1'; A=anymatrix(sA);
      case 510, sA='expm/fahi19r2'; A=anymatrix(sA);
      case 511, sA='expm/fahi19r3'; A=anymatrix(sA);
      case 512, sA='expm/fahi19r4'; A=anymatrix(sA);
      case 513, sA='expm/fasi7'; A=anymatrix(sA);
      case 514, sA='expm/jemc05r1'; A=anymatrix(sA);
      case 515, sA='expm/jemc05r2'; A=anymatrix(sA);
      case 516, sA='expm/kase99'; A=anymatrix(sA);
      case 517, sA='expm/kela89r1'; A=anymatrix(sA);
      case 518, sA='expm/kela89r2'; A=anymatrix(sA);
      case 519, sA='expm/kela98r1'; A=anymatrix(sA);
      case 520, sA='expm/kela98r2'; A=anymatrix(sA);
      case 521, sA='expm/kela98r3'; A=anymatrix(sA);
      case 522, sA='expm/kuda10'; A=anymatrix(sA);
      case 523, sA='expm/lara17r1'; A=anymatrix(sA);
      case 524, sA='expm/lara17r2'; A=anymatrix(sA);
      case 525, sA='expm/lara17r3'; A=anymatrix(sA);
      case 526, sA='expm/lara17r4'; A=anymatrix(sA);
      case 527, sA='expm/lara17r5'; A=anymatrix(sA);
      case 528, sA='expm/lara17r6'; A=anymatrix(sA);
      case 529, sA='expm/mopa03r1'; A=anymatrix(sA);
      case 530, sA='expm/mopa03r2'; A=anymatrix(sA);
      case 531, sA='expm/naha95'; A=anymatrix(sA);
      case 532, sA='expm/nies19'; A=anymatrix(sA);    
      case 533, sA='expm/pang85r1'; A=anymatrix(sA);
      case 534, sA='expm/pang85r2'; A=anymatrix(sA);
      case 535, sA='expm/pang85r3'; A=anymatrix(sA);
      case 536, sA='expm/ross8'; A=anymatrix(sA);
      case 537, sA='expm/trem05'; A=anymatrix(sA);
      case 538, sA='expm/tsin13'; A=anymatrix(sA);
      case 539, sA='expm/ward77r1'; A=anymatrix(sA);
      case 540, sA='expm/ward77r2'; A=anymatrix(sA); 
      case 541, sA='expm/ward77r3'; A=anymatrix(sA);
      case 542, sA='expm/ward77r4'; A=anymatrix(sA);      
    end   
end