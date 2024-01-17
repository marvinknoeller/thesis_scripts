function [X_in_between,ww,Pol_1,Pol_2,ts,der_X] = SetupFarField(X,M,num_x,mu_rel,eps_rel,ts)

if nargin == 6
    [~,~,coefs,~,ts] = splinepoints(X,M,ts);
else
    [~,~,coefs,~,ts] = splinepoints(X,M);
end
[X_in_between,der_X,derder_X,tt] = allpoints(coefs,ts,num_x,M);
delta_t = tt(2:end) - tt(1:end-1);
ww = zeros(1,size(X_in_between,2));
Pol_1 = zeros(3,3,size(X_in_between,2));
Pol_2 = zeros(3,3,size(X_in_between,2));
ww(1) = (delta_t(1)+delta_t(2))/6 * norm(der_X(:,1));
ww(end) = (delta_t(end)+delta_t(end-1))/6 *norm(der_X(:,end));
for ki = 2:2:(size(X_in_between,2)-1)
    ww(ki) = (delta_t(ki-1)+delta_t(ki))*4/6*norm(der_X(:,ki));   %Mid
    if ki<(size(X_in_between,2)-1)
        ww(ki+1) =  (delta_t(ki-1) + delta_t(ki) + delta_t(ki+1) + delta_t(ki+2))*1/6*norm(der_X(:,ki+1)); %End
    end
end
const_mu = 2/(mu_rel+1);
const_eps = 2/(eps_rel+1);
for ki = 1:size(X_in_between,2)
    tan = der_X(:,ki)/(norm(der_X(:,ki)));
    uno = [1;-1;0]; 
    %normal vector
    nor = uno - dotReal(uno,tan)*tan;
    nor = nor/norm(nor);  
    if norm(derder_X(:,ki))<1e-10
        nor = uno - dotReal(uno,tan)*tan;
        nor = nor/norm(nor);
    else
        ppp = derder_X(:,ki)/norm(derder_X(:,ki));
        nor = cross(cross(der_X(:,ki),ppp),der_X(:,ki));
        nor = nor/norm(nor);
    end
    %conormal vector
    bnor = cross(tan,nor);
    %the matrix V, we need to diag.
    V = [tan,nor,bnor];
    Pol_1(:,:,ki) = V*diag([1,const_mu,const_mu])*(V.');
    Pol_2(:,:,ki) = V*diag([1,const_eps,const_eps])*(V.');
end
end

