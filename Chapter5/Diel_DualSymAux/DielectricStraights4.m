clear all
close all
clc

lenvec = [16,22,28];
hvec = [4,6,8];
Nvec = [5,6,7];
for kk = 1 : 3
L = lenvec(kk);
Var.length = L;
N = Nvec(kk);

%-------------------------------------------------------
%% Variables for the scattering problem
Var.eps_rel = 50;
Var.mu_rel = 50;
Var.kappa = 1;
%% Variables for the optimization

aa = 1;
bb =  1.;
Var.N = N;
rho = .05;
Var.aa = aa * rho;
Var.bb = bb * rho;

Var.n = 10; %number of nodes
Var.M = 11; %number of points in each subsegment including first, last point
%% Initial guess of curve
rng(1234)
t = linspace(0,1,Var.n);
height = L;%lambda_point * 1000 * heightp *1e-9 * kappa2;%4.;
R = 1;
num_x = Var.n;
h = hvec(kk);
x0 = -R * ones(1,num_x) + 0.02*(2*rand(1,num_x)-1);
y0 = 0 * zeros(1,num_x) + 0.02*(2*rand(1,num_x)-1);
z0 = -(h*t -h/2);
p1 = [x0;y0;z0];

x0 = R * ones(1,num_x) + 0.02*(2*rand(1,num_x)-1);
y0 = 0 * zeros(1,num_x) + 0.02*(2*rand(1,num_x)-1);
z0 = (h*t -h/2);

p2 = [x0;y0;z0];

x0 = 0 * ones(1,num_x) + 0.02*(2*rand(1,num_x)-1);
y0 = 0 * zeros(1,num_x) + 0.02*(2*rand(1,num_x)-1);
z0 = (h*t -h/2);

p3 = [x0;y0;z0];

t = linspace(0.2,0.8,5);
x0 = R/2*cos(pi*t+pi)+ 0.02*(2*rand(1,5)-1) - R/2;
y0 = 0*zeros(1,5)+ 0.02*(2*rand(1,5)-1);
z0 = R*sin(pi*t+pi)-h/2;

bow = [x0;y0;z0];

t = linspace(0.2,0.8,5);
x0 = R/2*cos(pi*t+pi)+ 0.02*(2*rand(1,5)-1) + R/2;
y0 = 0*zeros(1,5)+ 0.02*(2*rand(1,5)-1);
z0 = -R*sin(pi*t+pi)+h/2;

bow2 = [x0;y0;z0];

pcurve = [p1,bow,p3,bow2,fliplr(p2)];
Var.n = size(pcurve,2);
%%
tt = linspace(0,1,length(pcurve));
alpha = 0*tt;
lambda3 = 0e-5;
lambda2 = 0.005;
lambda1 = 10;
Var.lambda1 = lambda1;
Var.lambda2 = lambda2;
Var.lambda3 = lambda3;
CurveOnlyDielectricBFGS(Var,pcurve,1,alpha,strcat('DielectricCurve4Len',num2str(L)));
end