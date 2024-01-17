function t_s_curvature = Curvature_Constraint(p,M,num_x,initial_length)
[~,~,coefs,~,ts] = splinepoints(p,M);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));
tangents = derder_p./((norm_p_d)) - dotReal(der_p,derder_p)./((norm_p_d).^3).*der_p;

integrand = 1./norm_p_d .* dotReal(tangents,tangents);
%%
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
%%
t_s_curvature = sum(ww.*integrand)/initial_length;
end

