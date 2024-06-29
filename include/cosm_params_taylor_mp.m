function [s, m] = cosm_params_taylor_mp(A, varargin)
%COSM_PARAMS_TAYLOR_MP  Getting the scaling parameter s and 
% approximant degree m from the multiprecision algorithm for the matrix
% cosine using Taylor approximants.
%
%   This code requires the Advanpix Multiprecision Computing Toolbox (see
%   www.advanpix.com).
%
%   [...] = cosm_params_taylor_mp(...,'precision',DIGITS) specifies the number of digits to be
%   used in the computation. Default is mp.Digits() if A is of class mp. The
%   computation is performed in single or double arithmetic if A is of class
%   single or double, respectively, and the precision is not specified.
%
%   [...] = cosm_params_taylor_mp(...,'epsilon',EPSILON) specifies the tolerance to be
%   used in evaluating the Taylor approximants. Default is machine epsilon 
%   of the precision of A if A is of class 'single' or 'double', mp.eps() 
%   if A is of class 'mp'.
%
%   [...] = cosm_params_taylor_mp(...,'maxdegree',MAXDEGREE) specifies the maximum degree of
%   the Taylor approximants. Default is 500.
%
%   Reference: 
%   Algorithm 4.4 in A. H. Al-Mohy, N. J. Higham, and X. Liu, Arbitrary 
%   Precision Algorithms for Computing the Matrix Cosine and its Fréchet 
%   Derivative. SIAM J. Matrix Anal. Appl., 43(1):233-256, 2022

%% Input Validation.
    p = inputParser;
    addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
    addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
    addParameter(p, 'maxdegree', 500, @(x)(x == round(x) && x > 0));
  
    parse(p,A,varargin{:});

    A = p.Results.A;
    digits = p.Results.precision;
    epsilon = p.Results.epsilon;
    maxdegree = p.Results.maxdegree;
    
    [n1, n] = size(A);
    if n1 ~= n, error('The matrix ``A'' must be square.'); end

    if nargout > 2
        error('This function returns at most two values.');
    end
    
%% Determine whether the computation will be single, double or mp.
    current_digits = mp.Digits();
    if ~isempty(digits) % working precision specified
        mp.Digits(digits);
        A = mp(A);
    else
        if isa(A, 'double')
            digits = 16;
        elseif isa(A, 'single')
            digits = 8;
        else
            digits = mp.Digits();
        end
        mp.Digits(digits);
    end
    if isempty(epsilon)
        epsilon = myeps(class(A));
    end
    curr_epsilon = epsilon;
           

    I = eye(n,'mp');
    B = A * A; % working with A^2;
    B_double = double(B);
    
    if ~any(B_double(:)) % if B is a zero matrix
        error('A^2 is a zero matrix.');
    end
    
    alphas = zeros(1,maxdegree);
    alphas(1) = norm(B_double, 1);
    
    degrees = opt_degrees(maxdegree); % vector of optimal degrees for current approximant.
    maxdegree = degrees(end);        % max degree
    
    i_deg = 1; % i_deg: index of current degree - we start from 1
    m = degrees(i_deg);           % current degree
    s = 0;
    Bsq_powers = zeros(n,n,ceil(sqrt(maxdegree))+1,'mp');
    Bsq_powers(:,:,1) = I;
    Bsq_powers(:,:,2) = B;
    Bsq_powers_lgth = 2; % current length of sequence B is 2 (B_0 = I, B_1 = B);
    Bsq_powers_double = double(Bsq_powers); % form double precision powers of B.
    
    
    abs_err = false; % use relative error bound
    if ~abs_err, compute_normcosm = true;  end
    
    % optional: compute the sum in t^{ch}_{2m} in the error bound using some extra precision
    ext_prec = true; % slight improves the accuracy, and brings no difference to speed. 
    prec_factor = 1.2; % the factor that determines the extra precision
    
    degree_found = false;
    delta_pre = Inf;

    % Prepare the factorials required in evaluating the bound - the 
    % initialization "factorials = factorial(2*mp(0:maxdegree));  
    % factorials_double = double(factorials);" is not used since its
    % execution time can be significant when maxdegree is large.
    % Instead, we initialize as following and compute more factorials 
    % only when needed.
    
    factorials = zeros(1,maxdegree+1,'mp');
    factorials_double = zeros(1,maxdegree+1); % powers of B in double precision for evaluating norm of cos(A/2^s)
    factorials(1:m+1) = factorial(2*mp(0:m));
    factorials_double(1:m+1) = double(factorials(1:m+1));
    
    tempnormcosm1 = Inf; % 'previous' approximation to the norm of cos(A/2^s) used in the function eval_bound
    [delta_nxt, Bsq_powers, tempnormcosm1] = eval_bound(Bsq_powers, alpha(s, m), m, ext_prec);
    if ~abs_err, curr_epsilon = epsilon * tempnormcosm1; end
    if delta_nxt < curr_epsilon, degree_found = true; end  % Awad (< instead of <=)
    
    order_k = 3;
    while ~degree_found && m < maxdegree 
         if abs(delta_pre) <= abs(delta_nxt)^order_k   % Awad
            s = s + 1;
            compute_normcosm = true; % s increased, recompute the 1-norm of cos(A/2^s)
        else
            i_deg = i_deg + 1; % m_i
            m_pre = m; % store the previous degree m
            m = degrees(i_deg);
            % update the factorials required in evaluating the bound
            factorials(m_pre+2:m+1) = factorial(2*mp(m_pre+1:m));
            factorials_double(m_pre+2:m+1) = double(factorials(m_pre+2:m+1));
        end
        delta_pre = delta_nxt;
        [delta_nxt, Bsq_powers, tempnormcosm1] = eval_bound(Bsq_powers, alpha(s, m), m, ext_prec);
        if ~abs_err, curr_epsilon = epsilon * tempnormcosm1; end
        if delta_nxt <= curr_epsilon, degree_found = true; end
    end
    
 %% Now s and m found, restore mp working precision.
    mp.Digits(current_digits);
      


%% -----------------------       SUBFUNCTIONS       -----------------------

    function [err_bnd, Bsq_powers, tempnormcosm] = eval_bound(Bsq_powers, alpha, m, ext_prec)
    %EVAL_BOUND   Error bound and current approx. of the norm of cos(A/2^s).
    %
    %   EVAL_BOUND(BSQ_POWERS, ALPHA, M, EXT_PREC) computes an upper bound 
    %   ERR_BND on the truncation error of the Taylor approximant of 
    %   degree M to COS(A/2^s) and an approximation TEMPNORMCOSM to its 
    %   norm using available powers of B stored in BSQ_POWERS. EXT_PREC 
    %   is a Boolean variable that determines whether we compute the 
    %   scalar coefficients in the error bound using some extra precision. 
    %   ALPHA represents the value alpha^*_m(B/4^s).
    
        nu = floor(sqrt(m));
        for i = (Bsq_powers_lgth+1):nu+1
            Bsq_powers(:,:,i) = Bsq_powers(:,:,i-1) * Bsq_powers(:,:,2);
        end
        % update the double precision powers
        Bsq_powers_double(:,:,Bsq_powers_lgth+1:nu+1) = double(Bsq_powers(:,:,Bsq_powers_lgth+1:nu+1));
        
        Bsq_pwlg_incrd = false; % Bsq_powers_lgth increased = false
        Bsq_pwlg_pre = Bsq_powers_lgth; 
        Bsq_powers_lgth = nu + 1;
        if Bsq_powers_lgth>Bsq_pwlg_pre, Bsq_pwlg_incrd = true; end
        
        digits_old = mp.Digits();
        if ext_prec, mp.Digits(prec_factor*digits_old); end
        alpha = sqrt(mp(alpha));
        t_ch =  sum(alpha.^(2*(0:m))./factorials(1:m+1)); 
        err_bnd = abs(cosh(alpha) - t_ch);
        mp.Digits(digits_old);
        if compute_normcosm % estimate and update 1-norm of cos(A/2^s) using potentially lower precision (double)
            reshaped_vector = reshape(4.^(-s*(0:Bsq_powers_lgth-1))./...
                factorials_double(1:Bsq_powers_lgth),[1,1,Bsq_powers_lgth]);
            approx = sum(bsxfun(@times,...
                                Bsq_powers_double(:,:,1:Bsq_powers_lgth),... % reshaped_matrix,...
                                reshaped_vector),3);
            tempnormcosm = norm(approx, 1);
            % for the same s turn off the estimation if the value of tempnormcosm is
            % converging (when Bsq_powers_lgth increased but diff is small)!
            if Bsq_pwlg_incrd && abs(tempnormcosm1 - tempnormcosm) / ...
                    abs(tempnormcosm) < 0.1
                compute_normcosm = false;
            end
        else
            tempnormcosm = tempnormcosm1;
        end
    end

    function [est,mv] = normBd(B,d)
    %NORMAM   Estimate of 1-norm of power of matrix.
    %   NORMBD(B,d) estimates norm(B^d,1).
    %
    %   If B has nonnegative elements the estimate is exact.
    %   [EST,MV] = NORMBD(B,d) returns the estimate EST and the number MV of
    %   matrix-vector products computed involving B or B^*.

    %   Reference: 
    %   [1] A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 
    %   31(3):970-989, 2009.
    %   [2] M. Fasi and N. J. Higham, A Arbitrary Precision Scaling and 
    %   Squaring Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 
    %   40(4):1233-1256, 2019
    
        nc = 1; % number of columns in the iteration matrix in NORMEST1
        n1 = size(B,1);
        if isequal(B,abs(B))
            e = ones(n1,1);
            for c=1:d         % for positive matrices only
                e = B'*e;
            end
            est = norm(e,inf);
            mv = d;
        else
            % Use powers of B in Bsq_powers in order to efficienlty compute
            % B^d * z.
            % mult_index stores indeces of powers of B for multiplications to be done later
            mult_index = zeros(1, Bsq_powers_lgth);  
            d_dec = d;
            lgth_dec = min(Bsq_powers_lgth, d+1); 
            while d_dec > 0
                mult_index(lgth_dec) = floor(d_dec / (lgth_dec-1));
                d_dec = mod(d_dec, lgth_dec-1);
                lgth_dec = min(lgth_dec - 1, d_dec+1);
            end
        [est,~,~,it] = normest1(@afun_power,nc);
        mv = it(2)*nc*d;
        end

        function prod = afun_power(flag,Z)
        %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

            if isequal(flag,'dim')
                prod = n1;
            elseif isequal(flag,'real')
                prod = isreal(B);
            else
                [n2,~] = size(Z);
                if n2 ~= n1
                    error('Dimension mismatch')
                end
                if isequal(flag,'notransp')
                    % do the matrix-vector multiplication according to mult_index
                    for i = find(mult_index(2:end))
                        for ii = 1:mult_index(i+1)  
                            Z = double(Bsq_powers(:,:,i+1)) * Z; % Bsq_powers starts from the identity
                        end
                    end
                    if mult_index(1) ~= 0
                        for ii = 1:mult_index(1)
                            Z = B * Z;
                        end
                    end
                elseif isequal(flag,'transp')
                    for i = find(mult_index(2:end))
                        for ii = 1:mult_index(i+1)
                            Z = double(Bsq_powers(:,:,i+1)') * Z;
                        end
                    end
                    if mult_index(1) ~= 0
                        for ii = 1:mult_index(1)
                            Z = B' * Z;
                        end
                    end
                end
                prod = Z;
            end
        end
    end

    function y = alpha(s, m)
    %ALPHA    Estimate 4^{-s}*max(||B^d||^(1/d), ||X^(d+1)||^(1/(d+1))).
    %   ALPHA(S,M) estimates 4^{-s}*max(||B^d||^(1/d), ||B^(d+1)||^(1/(d+1))) 
    %   in the 1-norm, where d is the largest integer such that d(d-1) < m+2.
        
        d_star = floor((1 + sqrt(4 * m + 5)) / 2);

        % Experimental strategy
        if alphas(d_star+1) == 0 % The second term is not known, we try the bound.
            if alphas(d_star) == 0
                alphas(d_star) = normBd(B_double, d_star)^(1/d_star);
            end
            known = find(alphas ~= 0);
            low = min(known);
            high = max(known);
            bin_counter = 0;
            found_upper_bound = false;
            while low < high
                if (low + high == d_star+1)
                    if (alphas(d_star) > alphas(low) * alphas(high))
                        found_upper_bound = true;
                        break;
                    end
                end
                if bin_counter
                    low = low + 1;
                else
                    high = high - 1;
                end
                bin_counter = mod(bin_counter + 1, 2);
            end
            if found_upper_bound
                %             fprintf('.\n');
                y = alphas(d_star) / 4^s;
                %             i = j;
                return
            else
                assert(alphas(d_star+1) == 0)
                alphas(d_star+1) = normBd(B_double, d_star+1)^(1/(d_star+1));
            end
        end
        [y,~] = max(alphas(d_star:d_star+1));
        y = y / 4^s;
    end

    function degs = opt_degrees(mmax)
    %OPT_DEGREES    Vector DEGS stores optimal degrees of Taylor
    %approximants to the cosine.
        degs = [1,    2,    4,    6,    9,   12,   16,   20,   25,...
                30,   36,   42,   49,   56,   64,   72,   81,   90,  100,...
                110,  121,  132,  144,  156,  169,  182,  196,  210,  225,...
                240,  256,  272,  289,  306,  324,  342,  361,  380,  400,...
                420,  441,  462,  484,  506,  529,  552,  576,  600,  625,...
                650,  676,  702,  729,  756,  784,  812,  841,  870,  900,...
                930,  961,  992, 1024, 1056, 1089, 1122, 1156, 1190, 1225,...
                1260, 1296, 1332, 1369, 1406, 1444, 1482, 1521, 1560, 1600,...
                1640, 1681, 1722, 1764, 1806, 1849, 1892, 1936, 1980, 2025,...
                2070, 2116, 2162, 2209, 2256, 2304, 2352, 2401, 2450, 2500];
        degs = degs(1:find(degs <= mmax, 1, 'last'));
    end  

 %% Multiprecision Computing Toolbox utility functions.  
    function e = myeps(curr_class)
    %COMP_EPS	Machine epsilon.
    %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
        if(strcmp(curr_class, 'mp'))
            e = mp('eps');
        else
            e = eps(curr_class)/2;
        end
    end

end