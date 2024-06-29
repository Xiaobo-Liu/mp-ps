function scal_err = scal_small_error(err, num_digs)
%SCAL_SMALL_ERROR   Rescals errors that are smaller than a tenth of the
%unit roundoff of the working precision.
% 
% NUM_DIGTS is the equivalent decimal digits of the working precision and
% ERR is the vector of errors before scaling.
% SCAL_ERR returns the vector of scaled errors.

u = 10^(-num_digs);
min_ee_digs = -log10(min(nonzeros(err)));

num = length(err);
scal_err = err;
for ii=1:num
    if err(ii)==0
        scal_err(ii) = u/10;
    else 
        if err(ii)<=u/10
            ee_digits = -log10(err(ii));
            scal_err(ii) = 10^(-num_digs - (ee_digits-(num_digs+1)) / ...
                (min_ee_digs-(num_digs+1)));
        end
    end       
end
end