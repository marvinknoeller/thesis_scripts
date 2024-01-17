function t_s_curvature = nPdashbP(p,M,num_x,tstuetz,initial_length,R,S,T)
[~,~,coefs,~,ts] = splinepoints(p,M,tstuetz);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));
nPdash = zeros(3,length(R));
nPdash(:,1) = (R(:,2) - R(:,1)) / (tt(2) - tt(1));
nPdash(:,end) = (R(:,end) - R(:,end-1)) / (tt(end) - tt(end-1));

for ki = 2:length(R)-1
    nPdash(:,ki) = (R(:,ki+1) - R(:,ki-1)) / (tt(ki+1) - tt(ki-1));
end

integrand = dotReal(nPdash,S).^2.*norm_p_d;
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