function  [Val,Length] = CurvePenalty(p,M,num_x,W)
[~,~,coefs,~,ts] = splinepoints(p,M);
[p_in_between,der_p,~,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));

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
Length = sum(ww.*norm_p_d);
if nargin ==4
    Val = (sum(ww.*norm_p_d)-W)^2;
else 
    Val = inf;
end
end