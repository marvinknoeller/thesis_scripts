function [points,t_spline,coef,br,t_stuetz] = splinepoints(p,M,t_stuetz)
% Input: p, the current parametrization
%        M, the number of points in between 2 points of p
%
% Output: splinepoints, all points on the spline
if nargin ==2
    t_stuetz = zeros(1,size(p,2));
    t_stuetz(1) = 0;
    for k = 1:size(p,2)-1
       t_stuetz(k+1) = norm(p(:,k+1) - p(:,k)) + t_stuetz(k);
    end
t_stuetz = t_stuetz./t_stuetz(end);
t_stuetz2 = zeros(1,size(p,2));
t_stuetz2(1) = 0;
t_stuetz2(end) = 1;
d = 0;
for kk = 1 : size(p,2) -1
    d = d + norm(p(:,kk+1) - p(:,kk));
end
for kk = 2:size(p,2)-1
    t_stuetz2(kk) = t_stuetz2(kk-1) + norm(p(:,kk) - p(:,kk-1))/d;
end
% t_stuetz = t_stuetz2;
norm(t_stuetz-t_stuetz2);
assert(norm(t_stuetz-t_stuetz2)<1e-12);
end
if size(p,1) == 3
y1 = p(1,:);
y2 = p(2,:);
y3 = p(3,:);
t_spline = [];
for k = 1:size(p,2)-1
    step = (t_stuetz(k+1) - t_stuetz(k))/(M-1);
    t_spline = [t_spline, t_stuetz(k) : step : t_stuetz(k+1)];  
    t_spline = t_spline(1,1:end-1);     %delete the last one. So there are no double points
end

t_spline(end+1) = t_stuetz(end);
%% This is for not a knot
f1 = spline(t_stuetz,y1,t_spline);
f2 = spline(t_stuetz,y2,t_spline);
f3 = spline(t_stuetz,y3,t_spline);

pp1 = spline(t_stuetz,y1);
pp2 = spline(t_stuetz,y2);
pp3 = spline(t_stuetz,y3);
%% This is for the natural spline
% pp1 = csape(t_stuetz,y1,'variational');
% pp2 = csape(t_stuetz,y2,'variational');
% pp3 = csape(t_stuetz,y3,'variational');
% 
% f1 = ppval(pp1,t_spline);
% f2 = ppval(pp2,t_spline);
% f3 = ppval(pp3,t_spline);
%%
coef(:,:,1) = pp1.coefs;
coef(:,:,2) = pp2.coefs;
coef(:,:,3) = pp3.coefs;
br(1,:) = pp1.breaks;
br(2,:) = pp2.breaks;
br(3,:) = pp3.breaks;
% points = [f1(t_spline); f2(t_spline); f3(t_spline)];
points = [f1; f2; f3];
elseif size(p,1)==1
    y1 = p(1,:);
    t_spline = [];
    for k = 1:size(p,2)-1
        step = (t_stuetz(k+1) - t_stuetz(k))/(M-1);
        t_spline = [t_spline, t_stuetz(k) : step : t_stuetz(k+1)];
        t_spline = t_spline(1,1:end-1);     %delete the last one. So there are no double points
    end
    
    t_spline(end+1) = t_stuetz(end);
    %% This is for not a knot
    f1 = spline(t_stuetz,y1,t_spline);
    
    pp1 = spline(t_stuetz,y1);
    %% This is for the natural spline
    % pp1 = csape(t_stuetz,y1,'variational');
    % pp2 = csape(t_stuetz,y2,'variational');
    % pp3 = csape(t_stuetz,y3,'variational');
    %
    % f1 = ppval(pp1,t_spline);
    % f2 = ppval(pp2,t_spline);
    % f3 = ppval(pp3,t_spline);
    %%
    coef(:,:,1) = pp1.coefs;
    br(1,:) = pp1.breaks;
    % points = [f1(t_spline); f2(t_spline); f3(t_spline)];
    points = [f1];
end