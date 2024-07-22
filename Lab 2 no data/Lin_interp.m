function [T_0, H, Temps_out] = Lin_interp(Temps,Lin_dist)
    coefs = polyfit(Lin_dist(2:end),Temps(2:end),1);
    Temps_out = polyval(coefs,Lin_dist);
    
    H = coefs(1);
    T_0 = Temps_out(1);
end

