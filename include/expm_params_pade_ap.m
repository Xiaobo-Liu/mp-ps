function [s, m] = expm_params_pade_ap(A, varargin)
%EXPM_PARAMS_PADE_AP  Getting the scaling parameter s and approximant
% degree m from the arbitrary precision algorithm for the matrix
% exponential using Pade approximants.
%
%   This code requires the Advanpix Multiprecision Computing Toolbox (see
%   www.advanpix.com).
%
%   [...] = expm_params_pade_ap(...,'precision',DIGITS) specifies the number of digits to be
%   used in the computation. Default is mp.Digits() if A is of class mp. The
%   computation is performed in single or double arithmetic if A is of class
%   single or double, respectively, and the precision is not specified.
%
%   [...] = expm_params_pade_ap(...,'epsilon',EPSILON) specifies the tolerance to be used to
%   evaluate the approximants. Default is machine epsilon of the precision
%   of A if A is of class 'single' or 'double', mp.eps() if A is of class 'mp'.
%
%   [...] = expm_params_pade_ap(...,'maxscaling',MAXSQRTM) specifies the maximum number of
%   matrix products allowed during the squaring phase. Default is 100.
%
%   [...] = expm_params_pade_ap(...,'maxdegree',MAXDEGREE) specifies the maximum degree of
%   the Pade approximants. Default is 500.
%
%   Reference:
%   Algorithm 4.1 in M. Fasi and N. J. Higham, An Arbitrary Precision Scaling
%   and Squaring Algorithm for the Matrix Exponentia. SIAM J. Matrix Anal.
%   Appl., 40(4):1233-1256, 2019

%% Parse and validate input.
    p = inputParser;
    addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0));
    addParameter(p, 'epsilon', [], @(x)(x > 0 && x < 1));
    addParameter(p, 'maxscaling', 100, @(x)(x == round(x) && x > 0));
    addParameter(p, 'maxdegree', 500, @(x)(x == round(x) && x > 0));
    addParameter(p, 'abserr', false, @(x)(x == true || x == false));

    parse(p,A,varargin{:});

    A = p.Results.A;
    digits = p.Results.precision;
    epsilon = p.Results.epsilon;
    maxscaling = p.Results.maxscaling;
    maxdegree = p.Results.maxdegree;
    abserr = p.Results.abserr;

    usetaylor = false;
    opt_deg_pade = @opt_degrees_diagonal;
    eval_pade_error = @scalar_error_diagonal;
    alpha_pade = @(s, m)(alpha(s, m, m));
    

    alphas = zeros(1,maxdegree);
    alphas(1) = norm(A,1);

    [n2, n] = size(A);
    if n ~= n2
        error('The matrix ``A'' must be square.');
    end

    %% Parse and validate output.
    if nargout > 2
        error('This function returns at most four values.');
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
    %     epsilon = epsilon;
    curr_epsilon = epsilon;
        
        X = A;
        X0 = X;

        % Create the vector of optimal degrees for current approximant.
        % Update maxdegree appropriately.
        degrees = opt_deg_pade(maxdegree);
        maxdegree = degrees(end);        % max degree
        currcost = 3;                    % current cost
        m = degrees(currcost);           % current degree

        %% Try all degrees lower than maxdegree.
        s = 0;

        Xsq_powers = zeros(n,n,ceil(sqrt(2*maxdegree))+1,class(X0));
        Xsq_powers(1:n+1:n^2) = 1;
        Xsq_powers(:,:,2) = X0 * X0;
        
        Xsq_powers_ll = 2;

        found_degree = false;
        p_factor = 1.1;
        tempnormexpm1 = Inf;
        compute_normexpm = true;
        cost_step = 1;

        ext_prec = true;
        [a, Xsq_powers, tempcondq, tempnormexpm1] =...
            eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
        if ~abserr
            curr_epsilon = epsilon * tempnormexpm1;
        end
        old_tempcondq = tempcondq;
        while(~isfinite(tempnormexpm1))
            disp(tempnormexpm1)
            currcost = currcost + cost_step;
            s = s + 1;
            compute_normexpm = true;
            m = degrees(currcost);
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
        end
        a_old = Inf;
        tempnormexpm1_old = 1;
        while ~found_degree && m < maxdegree && s < maxscaling
            normqinv_bound = epsilon^(-1/8);
            if a <= curr_epsilon && tempcondq < normqinv_bound
                found_degree = true;
            elseif tempcondq >= normqinv_bound ||...
                    (a_old > 1 && abs(a_old) < abs(a)^2)
                assert(~(tempcondq >= normqinv_bound && usetaylor))
                X = X / 2;
                s = s + 1;
                compute_normexpm = true;
            else
                currcost = currcost + cost_step;
                m = degrees(currcost);
            end
            a_old = a;
            tempnormexpm1_old = tempnormexpm1;
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, m), s, m, ext_prec);
            if ~abserr
                curr_epsilon = epsilon * tempnormexpm1;
            end
        end

        % Degree and mxm_ps are now fixed and we keep scaling until the bound
        % on the forward error is smaller than the unit roundoff.
        if ~found_degree
            [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
            if ~abserr
                curr_epsilon = epsilon * tempnormexpm1;
            end
            %         while (~isfinite(a) || a >= curr_epsilon || tempcondq >= normqinv_bound) && s < maxscaling
            while (~isfinite(a) || a >= curr_epsilon) && s < maxscaling
                X = X / 2;
                s = s + 1;
                compute_normexpm = true;
                [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, m, ext_prec);
                if ~abserr
                    curr_epsilon = epsilon * tempnormexpm1;
                end
            end
        end
     

    %% Restore mp working precision.
    mp.Digits(current_digits);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             SUBFUNCTIONS                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Scalar bound evaluation and Pade approximation.
    function [e, Xsq_powers, tempcondq, tempnormexpm] = scalar_error_diagonal(Xsq_powers, x, s, m, ext_prec)
    %SCALAR_ERROR_DIAGONAL   Error in Pade approx. to the exponential.
    %   SCALAR_ERROR_DIAGONAL(X,M) computes an upper bound on the
    %   truncation error of the [M/M] Pade approximant to EXP(A).
        digits_old = mp.Digits();
        if ext_prec
            mp.Digits(p_factor*digits_old);
        end
        m = mp(m);
        x = mp(x);
        l_num = mp(0:1:m);
        c_num = (factorial(m)/factorial(2 * m)) ./...
                (factorial(m-l_num).*factorial(l_num)) .* factorial(2*m-l_num);
        c_den = (-1).^l_num .* c_num;
        xx = cumprod([1, repmat(x, 1, m)]);

        [Ue, Xsq_powers] = polyvalm_ps(Xsq_powers, s, double(c_num(1:2:end)), 'double');
        if m >= 1
            [Uo, Xsq_powers] =...
                polyvalm_ps(Xsq_powers, s, double(c_num(2:2:end)), 'double');
            Uo = double(X0) / 2^s * Uo;
        else
            Uo = zeros(size(X0), 'like', Ue);
        end

        tempq = Ue - Uo;

        tempqinv = inv(double(tempq));
        tempnormqinv = norm(tempqinv, 1);
        tempcondq = norm(tempq, 1) * tempnormqinv;
        tempnormexpm = norm(2 * (tempq \ Uo) + eye(size(X0), 'like', tempq));

        e = tempnormqinv * abs(exp(x) * sum(c_den .* xx) - sum(c_num .* xx));
        mp.Digits(digits_old);

    end

    function [c,mv] = mynormAm(X,p)
    %NORMAM   Estimate of 1-norm of power of matrix.
    %   NORMAM(A,m) estimates norm(A^m,1).
    %   If A has nonnegative elements the estimate is exact.
    %   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
    %   matrix-vector products computed involving A or A^*.

    %   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
    %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
    %   970-989, 2009.

        tt = 1;
        n1 = length(X);
        if isequal(X,abs(X))
            e = ones(n1,1);
            for j=1:p         % for positive matrices only
                e = X'*e;
            end
            c = norm(e,inf);
            mv = p;
        else
            % Use powers of X in Xsq_powers in order to efficienlty compute
            % A^m * y.
            if usetaylor
                mult = zeros(1, Xsq_powers_ll);
                p_dec = p;
                ll_dec = min(Xsq_powers_ll, p+1);
                while p_dec > 0
                    mult(ll_dec) = floor(p_dec / (ll_dec-1));
                    p_dec = mod(p_dec,(ll_dec-1));
                    ll_dec = min(ll_dec - 1, p_dec + 1);
                end
            else
                ll = size(Xsq_powers, 3);
                mult = zeros(1, ll);
                p_dec = p;
                ll_dec = ll;
                while p_dec > 1 && ll_dec > 1
                    mult(ll_dec) = floor(p_dec / (2*(ll_dec-1)));
                    %     p_dec = p_dec - mult(ll_dec) * 2 * (ll_dec-1)
                    p_dec = mod(p_dec,2*(ll_dec-1));
                    ll_dec = min(ll_dec - 1, floor(p_dec/2) + 1);
                end
                mult(1) = p_dec;
            end

            [c,~,~,it] = normest1(@afun_power,tt);
            mv = it(2)*tt*p;
        end

        function Z = afun_power(flag,Y)
        %AFUN_POWER  Function to evaluate matrix products needed by NORMEST1.

            if isequal(flag,'dim')
                Z = n1;
            elseif isequal(flag,'real')
                Z = isreal(A);
            else
                [n2,~] = size(Y);
                if n2 ~= n1
                    error('Dimension mismatch')
                end
                if isequal(flag,'notransp')
                    for i = find(mult(2:end))
                        for ii = 1:mult(i+1)
                            Y = double(Xsq_powers(:,:,i+1)) * Y;
                        end
                    end
                    if mult(1) ~= 0
                        for ii = 1:mult(1)
                            Y = X * Y;
                        end
                    end
                elseif isequal(flag,'transp')
                    for i = find(mult(2:end))
                        for ii = 1:mult(i+1)
                            Y = double(Xsq_powers(:,:,i+1)') * Y;
                        end
                    end
                    if mult(1) ~= 0
                        for ii = 1:mult(1)
                            Y = X' * Y;
                        end
                    end
                end

                Z = Y;

            end

        end
    end

    function y = alpha(s, k, m)
    %ALPHA    Estimate max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1))).
    %   ALPHA(A,K,M) estimates max(||X^p||^(1/p), ||X^(p+1)||^(1/(p+1)))
    %   where p is the largest integer such that p(p-1) < K+M+1.
        j = floor((1 + sqrt(4 * (m + k) + 5)) / 2);

        % Experimental strategy
        if alphas(j+1) == 0 % The second term is not known, we try the bound.
            if alphas(j) == 0
                alphas(j) = mynormAm(double(X0), j)^(1/j);
            end
            known = find(alphas ~= 0);
            low = min(known);
            high = max(known);
            bin_counter = 0;
            found_upper_bound = false;
            while low < high
                if (low + high == j+1)
                    if (alphas(j) > alphas(low) * alphas(high))
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
                y = alphas(j) / 2^s;
                %             i = j;
                return
            else
                assert(alphas(j+1) == 0)
                alphas(j+1) = mynormAm(double(X0), j+1)^(1/(j+1));
            end
        end
        [y,i] = max(alphas(j:j+1));
        y = y / 2^s;
    end

    function [Y, Xsq_powers] = polyvalm_ps(Xsq_powers, s, c, outputclass)
    %POLYVALM_PS    Paterson-Stockmeyer method for matrix polynomials.
    %   POLYVALM_PS(A,C) evaluates the matrix polynomial
    %       Y = C(1)*EYE(N) + C(2)*A + ...  + C(M)*A^(M-1) + C(M+1)*A^M
    %   where A is an N-by-N matrix and C is a vector of length M+1
    %   using Paterson-Stockmeyer algorithm.

        if nargin == 3
            outputclass = class(X0);
        end
        mm = length(c) - 1;
        I = Xsq_powers(:,:,1);

        ss = ceil(sqrt(mm));
        if ss == 0
            Y = c(1) * Xsq_powers(:,:,1);
            return
        end
        rr = floor(mm/ss);

        % Compute first ss+1 powers.
        for i=Xsq_powers_ll+1:ss+1
            if mod(i,2) == 1 % even power
                i_half = (i - 1) / 2 + 1;
                Xsq_powers(:,:,i) = Xsq_powers(:,:,i_half)^2;
            else
                assert(all(all(Xsq_powers(:,:,i) == zeros(n,n))));
                Xsq_powers(:,:,i) = Xsq_powers(:,:,2) * Xsq_powers(:,:,i-1);
            end
        end
        Xsq_powers_ll = max(Xsq_powers_ll,ss+1);
       

        mpowers = Xsq_powers;
        for i = 1:ss+1
            mpowers(:,:,i) = mpowers(:,:,i)/(2^(2*s*(i-1)));
        end
        mpowers = cast(mpowers, outputclass);

        % Evaluate last polynomial of degree m mod ss.
        digits_old = mp.Digits();
        mp.Digits(1.2 * digits_old);
        B = mp(c(mm+1) * mpowers(:,:,mm-ss*rr+1));
        for j=mm-1:-1:ss*rr
            if j == ss*rr
                B = B + c(ss*rr+1)*I;
            else
                B = B + c(j+1) * mpowers(:,:,mm-ss*rr-(mm-j)+1);
            end
        end
        mp.Digits(digits_old);

        % Evaluate polynomials of degree ss-1 and evaluate main polynomial using
        % Horner's method.
        Y = cast(B, outputclass);
        for kk=rr-1:-1:0
            mp.Digits(1.2 * digits_old);
            B = zeros(size(X0,1), outputclass);
            B = B + c(ss*kk+1) * I;
            for j=1:ss-1
                B = B + c(ss*kk+j+1) * mpowers(:,:,j+1);
            end
            mp.Digits(digits_old);
            Y = Y * mpowers(:,:,ss+1) + cast(B, outputclass);
        end

    end

    function degs = opt_degrees_diagonal(nmax)
        degs = [1,    2,    3,    4,    6,    8,   10,   12,   15,...
                18,   21,   24,   28,   32,   36,   40,   45,   50,   55,...
                60,   66,   72,   78,   84,   91,   98,  105,  112,  120,...
                128,  136,  144,  153,  162,  171,  180,  190,  200,  210,...
                220,  231,  242,  253,  264,  276,  288,  300,  312,  325,...
                338,  351,  364,  378,  392,  406,  420,  435,  450,  465,...
                480,  496,  512,  528,  544,  561,  578,  595,  612,  630,...
                648,  666,  684,  703,  722,  741,  760,  780,  800,  820,...
                840,  861,  882,  903,  924,  946,  968,  990, 1012, 1035,...
                1058, 1081, 1104, 1128, 1152, 1176, 1200, 1225, 1250, 1275,...
                1300, 1326, 1352, 1378, 1404, 1431, 1458, 1485, 1512, 1540,...
                1568, 1596, 1624, 1653, 1682, 1711, 1740, 1770, 1800, 1830,...
                1860, 1891, 1922, 1953, 1984, 2016, 2048, 2080, 2112, 2145,...
                2178, 2211, 2244, 2278, 2312, 2346, 2380, 2415, 2450, 2485,...
                2520];
        degs = degs(1:find(degs >= nmax, 1, 'first'));
    end

    %% Multiprecision Computing Toolbox utility functions.

    function e = myeps(curr_class)
    %COMP_EPS    Machine epsilon.
    %   COMP_EPS(CLASS) computes the machine epsilon of class CLASS.
        if(strcmp(curr_class, 'mp'))
            e = mp('eps');
        else
            e = eps(curr_class);
        end
    end
end