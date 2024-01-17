function [ints] = Curv_Pen(p,M,num_x)
[~,~,coefs,~,ts] = splinepoints(p,M);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
norm_p_d = (sqrt(sum(der_p.^2,1)));
tangents = derder_p./((norm_p_d).^(3/2)) - dotReal(der_p,derder_p)./((norm_p_d).^(7/2)).*der_p;
delta_t = tt(2:end) - tt(1:end-1);
ww = zeros(1,size(p_in_between,2));
ww(1) = (delta_t(1)+delta_t(2))/6 ; 
assert(norm(delta_t(1) - delta_t(2))<1e-14);
ww(end) = (delta_t(end)+delta_t(end-1))/6; 
assert(norm(delta_t(end) - delta_t(end-1))<1e-14)

for j1 = 1 : num_x - 1 %for all subsegments
    seg_tt = tt((M-1)*(j1-1)+1 : (M-1)*j1+1);
    difs = seg_tt(2:end) - seg_tt(1:end-1);
    mean = 1/(M-1) * sum(difs); %M-1 subsubsegments
    assert(norm(difs-mean)<1e-14)
end


for ki = 2:2:(size(p_in_between,2)-1)
    ww(ki) = (delta_t(ki-1)+delta_t(ki))*4/6;   %Mid
    if ki<(size(p_in_between,2)-1)
        ww(ki+1) =  (delta_t(ki-1) + delta_t(ki) + delta_t(ki+1) + delta_t(ki+2))*1/6; %End
    end
end


%%
ints = sqrt(ww).*(tangents);
assert(length(tt) == (M-1)*(num_x-1)+1) 
% M points in each segment.
% But then, the control points are doubled.
% M is odd!
% For each line segment 
%int_t_j^{t_{j+1}} f(x) dx \approx \sum_{i=1}^M w_j f(x_j)
end

