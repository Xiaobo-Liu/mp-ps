function [s, m] = expm_params_taylor_ap(A, varargin)
%EXPM_PARAMS_TAYLOR_AP  Getting the scaling parameter s and approximant
% degree m from the arbitrary precision algorithm for the matrix
% exponential using Taylor approximants.
%
%   This code requires the Advanpix Multiprecision Computing Toolbox (see
%   www.advanpix.com).
%
%   [...] = expm_params_taylor_ap(...,'precision',DIGITS) specifies the number of digits to be
%   used in the computation. Default is mp.Digits() if A is of class mp. The
%   computation is performed in single or double arithmetic if A is of class
%   single or double, respectively, and the precision is not specified.
%
%   [...] = expm_params_taylor_ap(...,'epsilon',EPSILON) specifies the tolerance to be used to
%   evaluate the approximants. Default is machine epsilon of the precision
%   of A if A is of class 'single' or 'double', mp.eps() if A is of class 'mp'.
%
%   [...] = expm_params_taylor_ap(...,'maxscaling',MAXSQRTM) specifies the maximum number of
%   matrix products allowed during the squaring phase. Default is 100.
%
%   [...] = expm_params_taylor_ap(...,'maxdegree',MAXDEGREE) specifies the maximum degree of
%   the Taylor approximants. Default is 500.
%
%   Reference:
%   Algorithm 4.1 in M. Fasi and N. J. Higham, An Arbitrary Precision Scaling
%   and Squaring Algorithm for the Matrix Exponentia. SIAM J. Matrix Anal.
%   Appl., 40(4):1233-1256, 2019

%% Parse and validate input.
    p = inputParser;
    addRequired(p, 'A', @(x)(ismatrix(x) && isnumeric(x)));
    addParameter(p, 'precision', [], @(x)(x == round(x) && x > 0))
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

    usetaylor = true;
    opt_deg_pade = @opt_degrees_taylor;
    eval_pade_error = @scalar_error_taylor;
    alpha_pade = @(s, m)(alpha(s, m, 0));
  
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

        factorials = factorial(mp(0:maxdegree));
        factorials_double = double(factorials);
        Xsq_powers = zeros(n,n,ceil(sqrt(maxdegree))+1,'double');
        Xsq_powers(1:n+1:n^2) = 1;
        Xsq_powers(:,:,2) = X0;
      
        Xsq_powers_ll = 2;

        found_degree = false;
        p_factor = 1.1;
        tempnormexpm1 = Inf;
        compute_normexpm = true;
        cost_step = 1;

        if norm(A,1) > 1e7
            ext_prec = false;
            [a, Xsq_powers, tempcondq, tempnormexpm1] =...
                eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
            while (a > epsilon * tempnormexpm1 || ~isfinite(tempnormexpm1)) && s <= maxscaling
                s = s + 1;
                X = X / 2;
                compute_normexpm = true;
                [a, Xsq_powers, tempcondq, tempnormexpm1] = eval_pade_error(Xsq_powers, alpha_pade(s, maxdegree), s, maxdegree, ext_prec);
            end
        end

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
            bound_time = tic;
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

    function [e, Xsq_powers, tempcondq, tempnormexpm] = scalar_error_taylor(Xsq_powers, x, s, m, ext_prec)
    %SCALAR_ERROR_TAYLOR   Error in Taylor approx. to the exponential.
    %   SCALAR_ERROR_TAYLOR(A,X,M) computes an upper bound on the
    %   truncation error of the Taylor approximant of degree M to EXP(A).

        ss = ceil(sqrt(m));
        for i = Xsq_powers_ll+1:ss+1
            Xsq_powers(:,:,i) = Xsq_powers(:,:,i-1) * Xsq_powers(:,:,2);
        end
        Xsq_powers_ll = ss+1;
   

        digits_old = mp.Digits();
        if ext_prec
            mp.Digits(p_factor*digits_old);
        end
        x = mp(x);
        y = sum(x.^(0:m)./factorials(1:m+1));
        e = abs(y - exp(x));
        mp.Digits(digits_old);
        tempcondq = 1;
        if compute_normexpm
            reshaped_vector = reshape(2.^(-s*(0:Xsq_powers_ll-1))./factorials_double(1:Xsq_powers_ll),[1,1,Xsq_powers_ll]);
            approx = sum(bsxfun(@times,...
                                Xsq_powers(:,:,1:Xsq_powers_ll),... % reshaped_matrix,...
                                reshaped_vector),3);
            tempnormexpm = norm(approx, 1);
            if abs(tempnormexpm1 - tempnormexpm) / abs(tempnormexpm) < sqrt(eps)
                compute_normexpm = false;
                %             disp('.');
            end
        else
            tempnormexpm = tempnormexpm1;
        end
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
            mult = zeros(1, Xsq_powers_ll);
            p_dec = p;
            ll_dec = min(Xsq_powers_ll, p+1);
            while p_dec > 0
                mult(ll_dec) = floor(p_dec / (ll_dec-1));
                p_dec = mod(p_dec,(ll_dec-1));
                ll_dec = min(ll_dec - 1, p_dec + 1);
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
  
    function degs = opt_degrees_taylor(nmax)
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
        degs = degs(1:find(degs <= nmax, 1, 'last'));
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