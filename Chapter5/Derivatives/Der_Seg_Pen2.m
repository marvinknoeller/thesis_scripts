function [der_seg] = Der_Seg_Pen2(p,h,M,num_x,initial_Len)

[~,~,coefs,~,ts] = splinepoints(p,M);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));
[~,~,coefsh,~,tsh] = splinepoints(h,M,ts);
[~,der_h,derder_h,~] = allpoints(coefsh,tsh,num_x,M);

for ite = 1:num_x-1 %loop over all segments
    points_seg = tt(M*(ite-1)+1-(ite-1):M*ite-(ite)+1);
    dp_seg = der_p(:,M*(ite-1)+1-(ite-1):M*ite-(ite)+1);
    dh_seg = der_h(:,M*(ite-1)+1-(ite-1):M*ite-(ite)+1);
    
    norm_p_seg = (sqrt(sum(dp_seg.^2,1)));
    delta_t_seg = points_seg(2:end)-points_seg(1:end-1);
    ww_seg = zeros(1,size(points_seg,2));
    ww_seg(1) = (delta_t_seg(1)+delta_t_seg(2))/6 ;
    ww_seg(end) = (delta_t_seg(end)+delta_t_seg(end-1))/6 ;
    for ki = 2:2:(size(points_seg,2)-1)
        ww_seg(ki) = (delta_t_seg(ki-1)+delta_t_seg(ki))*4/6;   %Mid
        if ki<(size(points_seg,2)-1)
            ww_seg(ki+1) =  (delta_t_seg(ki-1) + delta_t_seg(ki) + delta_t_seg(ki+1) + delta_t_seg(ki+2))*1/6; %End
        end
    end
    int_seg1(ite) = sum(ww_seg.*norm_p_seg);
    int_seg2(ite) = sum(ww_seg.*(dotReal(dp_seg,dh_seg)./(norm_p_seg)));
    int_seg_der3(ite) = 2*int_seg1(ite)*int_seg2(ite);
end

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
der_seg = sum(-2*int_seg2 .* (1/(num_x-1)*initial_Len - int_seg1)) / (initial_Len^2) ;
end