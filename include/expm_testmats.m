function [A, n_mats] = expm_testmats(k, n)
%EXPM_TESTMATS Test matrices for the matrix exponential that are divided
%into six groups.
%
%   [A, NMATS] = EXPM_TESTMATS(K,N) selects the K'th test matrix.
%   NMATS is the number of test matrices available.
%   N sets the dimension of the matrices if the size is undetermined.
%
%   In total, there are 30 + 17 + 28 + 2 + 40 + 20  = 137 test matrices.
%   Matrices of determined dimension have size ranging from 2-by-2 to 
%   21-by-21, e.g., for the default value n = 20 (if n is unspecified), 
%   the testset contain square matrices of size 2 to 21.
%
%   Groups 1-4 consists of matrices that are more general and mostly from 
%   the MATLAB gallery, and they are classified by having
%   i) positive real eigenvalues, ii) real eigenvalues, 
%   iii) complex eigenvalues, and iv) zero as an eigenvalue.
% 
%   Group 5-6 are special variants of the matrices in Group 1-4 or 
%   thoughtfully constructed matrices that have been used in the literature.

    indices = [101:130, 201:217, 301:328, 401:402, 501:540, 601:620];

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

        %% Group 1: positive real eigenvalues
      case 101, sA='cauchy';A=gallery(sA,n); %symm % optional parameters
      case 102, sA='condex';A=gallery(sA,2,4,6);
      case 103, sA='condex';A=gallery(sA,3,2);
      case 104, sA='condex';A=gallery(sA,n,3); %lower triangular
      case 105, sA='condex';A=gallery(sA,n,4,100); %symm real
      case 106, sA='dorr';A=full(gallery(sA,n,100));
      case 107, sA='dramadah';A=gallery(sA,n,2);
      case 108, sA='frank';A=gallery(sA,n);
      case 109, sA='gcdmat';A=gallery(sA,n); %symm real
      case 110, sA='grcar';A=gallery(sA,n);

      case 111, sA='hanowa';A=gallery(sA,n);
      case 112, A=hilb(n); %symm real
      case 113, sA='invhess';A=gallery(sA,n);
      case 114, sA='jordbloc';A=gallery(sA,n,1); %symm real
      case 115, sA='kahan';A=gallery(sA,n);
      case 116, sA='lehmer';A=gallery(sA,n); %symm real
      case 117, sA='minij';A=gallery(sA,n); %symm real
      case 118, sA='moler';A=gallery(sA,n); %symm real
      case 119, sA='parter';A=gallery(sA,n);
      case 120, sA='pei';A=gallery(sA,n);

      case 121, sA='poisson';A=gallery(sA,ceil(sqrt(n))); %symm real, n^2
      case 122, sA='prolate';A=gallery(sA,n,1); %symm real Toeplitz
      case 123, sA='randcorr';A=gallery(sA,n); %symm real
      case 124, sA='sampling';A=gallery(sA,n);
      case 125, sA='toeppd';A=gallery(sA,n); %symm real
      case 126, sA='tridiag';A=full(gallery(sA,n));
      case 127, sA='wathen';A=full(gallery(sA,ceil(n^(1/4)),ceil(n^(1/4)))); %symm real
      case 128, sA='wilk';A=full(gallery(sA,3));
      case 129, sA='wilk';A=full(gallery(sA,4));
      case 130, sA='wilk';A=full(gallery(sA,5)); %symm real

        %% Group 2: real eigenvalues.
      case 201, sA='binomial'; A=gallery(sA,n);
      case 202, sA='fiedler';A=gallery(sA,n);
      case 203, sA = 'house'; [v, beta] = gallery(sA, n);
        A = eye(n) - beta * (v * v');
      case 204, sA='jordbloc';A=gallery(sA,n,2); %symm real
      case 205, sA='kms';A=gallery(sA,n); %symm real
      case 206, sA='lesp';A=gallery(sA,n);
      case 207, sA='lotkin';A=gallery(sA,n);
      case 208, sA='orthog';A=gallery(sA,n,1); %symm real
      case 209, sA='orthog';A=gallery(sA,n,2); %symm real
      case 210, sA='orthog';A=gallery(sA,n,5); %symm real

      case 211, sA='orthog';A=gallery(sA,n,6); %symm real
      case 212, sA='orthog';A=gallery(sA,n,-1); %symm real
      case 213, sA='redheff';A=gallery(sA,n);
      case 214, sA='riemann';A=gallery(sA,n);
      case 215, sA='ris';A=gallery(sA,n,1e1);
      case 216, sA='wilk';A=full(gallery(sA,21)); %symm real
      case 217, sA='clement';A=gallery(sA,n,1); %symm


        %% Group 3: complex eigenvalues
      case 301, sA='chebspec';A=gallery(sA,n);
      case 302, sA='chebvand';A=gallery(sA,n);
      case 303, sA='chow';A=gallery(sA,n);
      case 304, sA='circul';A=gallery(sA,n);
      case 305, sA='cycol';A=gallery(sA,n);
      case 306, sA='dramadah';A=gallery(sA,n,1);
      case 307, sA='dramadah';A=gallery(sA,n,3);
      case 308, sA='forsythe';A=gallery(sA,n);
        %   case 309, sA='invol';A=gallery(sA,n); % causing overflow with n=500, removed
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

        %% Group 4: contain zero as an eigenvalue
      case 401, sA='gearmat';A=gallery(sA,n);
      case 402, sA='neumann';A=gallery(sA,ceil(sqrt(n))^2);
      
        %% Group 5: from the literature, used by Awad and Nick in expm.
      case 501, A = [4 2 0; 1 4 1; 1 1 4]; % \cite[Test 1]{ward77}
      case 502, A = [29.87942128909879 0.7815750847907159 -2.289519314033932
             .7815750847907159   25.72656945571064    8.680737820540137
             -2.289519314033932   8.680737820540137  34.39400925519054]; % \cite[Test 2]{ward77}
      case 503, A = [-131 19 18; -390 56 54; -387 57 52]; % \cite[Test 3]{ward77}
      case 504, A = gallery('forsythe',10,1e-10,0); % \cite[Test 4]{ward77}
      case 505, T = [1 10 100; 1 9 100; 1 11 99]; 
                A = T*[0.001 0 0; 0 1 0; 0 0 100] / T; % \cite[p. 370]{naha95}
      case 506, A = [0.1 1e6; 0 0.1]; % \cite[Ex.~2]{kela98}.
      case 507, A = [0  3.8e3 0    0   0
                     0 -3.8e3 1    0   0
                     0 0     -1  5.5e6 0
                     0 0      0 -5.5e6 2.7e7
                     0 0      0   0   -2.7e7]; % \cite[p.~655]{kela98}
      case 508, w = 1.3; x = 1e6; m = 8;
                A = (1/m) * [w*ones(m/2) x*ones(m/2)
                     zeros(m/2)  -w*ones(m/2)]; % \cite[Ex.~3.10]{dipa00}
      case 509, A = rosser; A = 2.05*A/norm(A,1);  % Bad case for expm re. cost
      case 510, A = [0 1e4; -1e4 0];  % exp = [cos(x) sin(x); - sin(x) cos(x)], x = 100;
      case 511, A = 1e2*triu(randn(n),1);  % Nilpotent.
      case 512, A = zeros(n); A(n+1:n+1:n^2) = 1:n-1; % log of Cholesky factor of Pascal matrix. \cite{edst03}.
      case 513, A = [48 -49 50 49; 0 -2 100 0; 0 -1 -2 1; -50 50 50 -52]; % \cite[p.~206]{kela89}
      case 514, A = [0    30 1   1  1  1
                    -100   0 1   1  1  1
                     0     0 0  -6  1  1
                     0     0 500 0  1  1
                     0     0 0   0  0  200
                     0     0 0   0 -15 0]; % \cite[p.~7, Ex I]{pang85}
      case 515, A = gallery('triw',n,1);  m = (n-1)/2;
                A = A - diag(diag(A)) + diag(-m:m)*sqrt(-1);
                for i = 1:n-1, A(i,i+1) = -2*(n-1)-2 + 4*i; end 
                % Max Fasi's interpretation of their matrix for arbitrary n.
                % n = 31 corresponds to the matrix in \cite[p.~9, Ex II]{pang85}    
      case 516, A = gallery('triw',n,1,1);
                A = A - diag(diag(A)) + diag(-(n-1)/2:(n-1)/2); % \cite[p.~10, Ex III]{pang85}

      case 517, A = [0 1e6; 0 0];   % Same as case 506 but with ei'val 0.1 -> 0 % \cite[Ex.~5]{kela89}
      case 518, g = [0.6 0.6 4.0]; b = [2.0 0.75];
                 A = [-g(1)       0    g(1)*b(1)
                        0        -g(2) g(2)*b(2)
                      -g(1)*g(3)  g(3) -g(3)*(1-g(1)*b(1))]; % \cite[(52)]{jemc05}
      case 519, g = [1.5 0.5 3.0 2.0 0.4 0.03]; b = [0.6 7.0];
                A1 = [-g(5)     0      0
                        0      -g(1)    0
                        g(4)     g(4)   -g(3)];
                A2 = [-g(6)       0    g(6)*b(2)
                        0        -g(2)  g(2)*b(1)
                        0         g(4) -g(4)];
                A = [zeros(3) eye(3); A2 A1]; % \cite[(55)]{jemc05}
      case 520, A = [-1 1e7; 0 -1e7]; % \cite[Ex.~3]{kela98}
      case 521, Thalf = [3.8235*60*24 3.10 26.8 19.9]/60;  % Half lives in seconds
                a = log(2)./Thalf;  % decay constant
                A = diag(-a) + diag(a(1:end-1),-1); % \cite[(21)]{mopa03}
      case 522, a1 = 0.01145; a2 = 0.2270;
                A = [-a1              0  0
                      0.3594*a1     -a2  0
                      0.6406*a1      a2  0]; % \cite[(26)]{mopa03}
      case 523, a = [4.916e-18
                     3.329e-7
                     8.983e-14
                     2.852e-13
                     1.373e-11
                     2.098e-6
                     9.850e-10
                     1.601e-6
                     5.796e-8
                     0.000];
                A = diag(-a) + diag(a(1:end-1),-1); % \cite[Table 1]{kase99}
      case 524, lambda = 1e6 * 1i; mu = 1/2*(-1+sqrt(1+4*lambda));
                A = [ 0, 1; lambda, -1 ] - mu*eye(2); % by Jitse Niesen
      case 525, A = [1 1e17;0 1]; % by Awad
      case 526, b = 1e3; x = 1e10;
                A = [ 1-b/2   b/2;       -b/2       1+b/2];
                A = [ A       x*ones(2);  zeros(2) -A]; % by Awad    
      case 527, b = 1e4; A = [1-b/2   b/2 ; -b/2   1+b/2]; % by Awad         
      case 528, b = 1e2; A = [1-b/2   b/2 ; -b^4/2   1+b/2]; % by Awad
      case 529, godunov_demo = [ 289  2064   336   128    80     32    16
                                1152    30  1312   512   288    128    32
                                 -29 -2000   756   384  1008    224    48
                                 512   128   640     0   640    512   128
                                1053  2256  -504  -384  -756    800   208
                                -287   -16  1712  -128  1968    -30  2032
                               -2176  -287 -1565  -512  -541  -1152  -289 ];
                A = godunov_demo/100; % from EIGTOOL
      case 530, A = 10*[0 1 2; -0.01 0 3; 0 0 0]; % \cite[(14.17), p. 141]{trem05}
      case 531, A = triu(schur(gallery('invol',13),'complex'),1);
      case 532, alpha = 1; beta = 1;  % No values are given in the paper % \cite{kuda10}
                A = -eye(n) + alpha/2*(diag(ones(n-1,1),1) + diag(ones(n-1,1),-1));
                A(1,2) = beta; A(n,n-1) = beta;
      case 533, A = [-3.328853448977761e-07 4.915959875924379e-18;
                      0    -4.915959875924379e-18]; % \cite[Benchmark #1]{lara17} \cite[Problem 1]{zhao17}
      case 534, A = [-2.974063693062615e-07      0   1.024464026382002e-14;
                      2.974063693062615e-07 -1.379680196333551e-13    0;
                      0      0    -1.024464026382002e-14]; % \cite[Benchmark #2]{lara17} \cite[Problem 2]{zhao17}
      case 535, A = [-2.421897905520424e-03            0      5.383443102348909e-03;
                      0       -3.200125487349701e-04            0; 
                      0        3.200125487349701e-04 -5.398342527725431e-03]; % \cite[Benchmark #3]{lara17} \cite[Problem 3]{zhao17}
      case 536, A = [-1.000000000000312e-04         0                  0  0;
                      1.000000000000000e-04 -1.000000000009379e-04     0  0;
                      0  1.000000000000000e-04 -1.188523972153541e-06  0;
                      0     0   1.188523972153541e-06 -1.024464026382002e-14]; % \cite[Benchmark #4]{lara17} \cite[Problem 4]{zhao17}
      case 537, A = sparse(...
                    [ 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6,...
                      7, 7, 8, 8, 9, 9,10,10,11,11,12,12],...
                    [ 1, 7, 1, 2, 2, 3,10, 3, 4, 4, 5,12, 5, 6, 6,...
                      7, 7, 8, 6, 9, 8,10,10,11,11,12],...
                    [ -1.000000000000049e-04
                       5.880666420493406e-14
                       1.000000000000000e-04
                      -4.926419193745169e-04
                       4.926419193745169e-04
                      -3.405151448232769e-06
                       2.980258985838552e-12
                       3.405151448232769e-06
                      -1.000000009110124e-04
                       1.000000000000000e-04
                      -1.000000033477380e-04
                       1.212838692746004e-09
                       1.000000000000000e-04
                      -1.000015370544945e-04
                       1.000000000000000e-04
                      -1.000000000588073e-04
                       1.000000000000000e-04
                      -3.885005720114481e-05
                       1.537023753355886e-09
                      -5.077325179294990e-11
                       3.885005720114481e-05
                      -1.000000029802590e-04
                       1.000000000000000e-04
                      -1.906345381077957e-05
                       1.906345381077957e-05
                      -1.212838692746004e-09]); A = full(A); % \cite[Benchmark #5]{lara17} \cite[Problem 5]{zhao17}
      case 538, A = sparse([1  1  2  2  2  3  3  3  4  5  5  6  6  7  7  8  8],...
                    [1  4  1  2  4  2  3  4  4  4  5  5  6  6  7  7  8],...
                    [ -2.930607054625170e-05
                       1.292290622141271e-07
                       2.446793135977101e-05
                      -2.106574217602557e-05
                       2.051479647948103e-08
                       2.106574217602557e-05
                      -9.549786402447881e-15
                       1.855074206409039e-12
                      -1.100000000000049e-04
                       1.000000000000000e-04
                      -4.926419193745169e-04
                       4.926419193745169e-04
                      -3.405151448232769e-06
                       3.405151448232769e-06
                      -1.000000091101239e-05
                       1.000000000000000e-05
                      -3.347737955438215e-12]); A = full(A); % \cite[Benchmark #6]{lara17} \cite[Problem 6]{zhao17}
       case 539, e = -1; f = -1.1;
                 A = blkdiag(compan([1, -4*e, 6*e^2, -4*e^3 e^4]), ...
                             compan([1, -3*f, 3*f^2, -f^3])); % 7-by-7 easy example. Idea from \cite[Matrix F, p. 81]{kags77a}
       case 540, A = gallery('triw',4,2^(60)) - diag([17 17 2 2]); % very ill conditioned \cite[Ex.~3]{dahi03}
                 A = A/1e4;  % Make less ill conditioned.

        %% Group 6: designed by Max Fasi, used by Max and Nick in mpexpm.
       case 601, A = full(gallery('tridiag',2,-1,-2,1));
       case 602, A = [1 1; 1 1 + 10 * eps];
       case 603, A = [10 0 0; 0 1 1; 0 1 1 + 10 * eps];
       case 604, A = []; for i = 1:4, A = blkdiag(A,...
                        full(gallery('tridiag',i,0,i,1))); end
       case 605, m = 10; D = diag([zeros(1, m-1) 1]); Q = triu(ones(m)); A = Q * D / Q;
       case 606, A = diag([0, 1, 1e6]);
       case 607, A = [1e-4 0; 0 1e4];
       case 608, A = ones(2,2);
       case 609, A = toeplitz([16-3i, (4+3i)/8, zeros(1,n-2)],[16-3i, -5, zeros(1,n-2)]);
       case 610, A = [1 1; 0 1e2];
       case 611, A = [1 1e3; 1e3 1];
       case 612, A = [1 2 3; 1 2 3; 1 2 3];
       case 613, t = -pi/2; A = [cos(t) -sin(t); sin(t) cos(t)];
       case 614, v = eye(n,1); A = eye(n) - v*v';
       case 615, A = [100 2 3; 4 5 6; 7 8 100];
       case 616, A = [1 1 1; 1 1 1 + 10*eps; 1 1 1+100*eps()];
       case 617, A = [1 2 3; 4 5 6; 7 8 1e2];
       case 618, A = [1 1 1 0.1; 1 1 1 10*eps; 1 1 1 100*eps; 1 1 1 1000*eps];
       case 619, A = gallery('randsvd', 8, 1e14);  
       case 620, A = gallery('randjorth',4,4,1e8) / 1e3;
    end   
end