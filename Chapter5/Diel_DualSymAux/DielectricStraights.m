clear all
close all
clc

lenvec = [4,6,8];
Nvec = [4,5,6];
% frVec = [400 500 600 700];
for kk = 1 : 3
    L = lenvec(kk);
    Var.length = L;
    N = Nvec(kk);

%-------------------------------------------------------
%% Variables for the scattering problem
Var.eps_rel = 5;%lambda_eps;
Var.mu_rel = 1;
Var.kappa = 1;
%% Variables for the optimization

aa = 1;
bb =  1.;
Var.N = N;
rho = .05;
Var.aa = aa * rho;
Var.bb = bb * rho;

Var.n = 20; %number of nodes
Var.M = 11; %number of points in each subsegment including first, last point
%% Initial guess of curve
rng(1234)
tt = linspace(0,1,Var.n);
% kappa2 = 2*pi/(lambda_point*1e-6);
height = L;%lambda_point * 1000 * heightp *1e-9 * kappa2;%4.;
radius = 0.02;%lambda_point* 1000 * radiusp *1e-9 * kappa2;
x0 =  2*(rand(1,length(tt)) - 1/2) * radius;
y0 =  2*(rand(1,length(tt)) - 1/2) * radius;
z0 = height* tt - height/2;

pcurve = [x0;y0;z0];
%%
tt = linspace(0,1,Var.n);
alpha = 0*tt;
lambda3 = 0e-5;
lambda2 = 0.0005;
lambda1 = 5;
Var.lambda1 = lambda1;
Var.lambda2 = lambda2;
Var.lambda3 = lambda3;
CurveOnlyDielectricBFGS(Var,pcurve,1,alpha,strcat('DielectricCurve1Len',num2str(L)));
end