function eescal = scal_small_error(ee_old, num_digs)
% EESCAL rescals too small errors to be no smaller than 10^{-NUM_DIGS-1} 
% in the vector EE 

u = 10^(-num_digs);
min_ee_digs = -log10(min(nonzeros(ee_old)));

num = length(ee_old);
eescal = ee_old;
for ii=1:num
    if ee_old(ii)==0
        eescal(ii) = u/10;
    else 
        if ee_old(ii)<=u/10
            ee_digits = -log10(ee_old(ii));
            eescal(ii) = 10^(-num_digs - (ee_digits-(num_digs+1)) / (min_ee_digs-(num_digs+1)));
        end
    end       
end
end