function [der_curve] = Der_Curvature_Constraint(p,h,M,num_x,initial_length)
[~,~,coefs,~,ts] = splinepoints(p,M);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));
[~,~,coefsh,~,tsh] = splinepoints(h,M,ts);
[~,der_h,derder_h,~] = allpoints(coefsh,tsh,num_x,M);

 integrand_der = 2*dotReal(derder_p,derder_h)./(norm_p_d.^3) ... 
     -3*dotReal(der_p,der_h)./(norm_p_d.^5) .* dotReal(derder_p,derder_p) ...
     -4*dotReal(der_p,derder_p)./(norm_p_d.^5) .* (dotReal(der_p,derder_h) + dotReal(der_h,derder_p)) ...
     +10*dotReal(der_p,der_h)./(norm_p_d.^7) .*dotReal(der_p,derder_p).^2 ...
     +2*dotReal(der_p,derder_p)./(norm_p_d.^5).*(dotReal(der_h,derder_p) + dotReal(der_p,derder_h)) ...
     -5./(norm_p_d.^7).*dotReal(der_p,derder_p).^2 .*dotReal(der_p,der_h);
 
%% weights
delta_t = tt(2:end) - tt(1:end-1);
ww = zeros(1,size(p_in_between,2));
    ww(1) = (delta_t(1)+delta_t(2))/6 ;
    ww(end) = (delta_t(end)+delta_t(end-1))/6 ;
    for ki = 2:2:(size(p_in_between,2)-1)
        ww(ki) = (delta_t(ki-1)+delta_t(ki))*4/6;   %Mid
        if ki<(size(p_in_between,2)-1)
            ww(ki+1) =  (delta_t(ki-1) + delta_t(ki) + delta_t(ki+1) + delta_t(ki+2))*1/6; %End
        end
    end

der_curve = sum(ww.*(integrand_der)) / initial_length;
end

