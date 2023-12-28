function coef = comput_coef_Pade(m, num_digs, options)
% compute the vector of Pade coefficients for the matrix exponential.

mp.Digits(num_digs);
m = mp(m);
l_num = mp(0:1:m);

if isequal(options,'numerator')
    coef = (factorial(m)/factorial(2 * m)) ./...
                (factorial(m-l_num).*factorial(l_num)) .*...
                factorial(2*m - l_num);
end

if isequal(options,'denominator')
    coef = (-1).^l_num .* ...
        (factorial(m)/factorial(2 * m)) ./ ...
                (factorial(m-l_num).*factorial(l_num)) .* ...
                factorial(2*m - l_num);
end
end