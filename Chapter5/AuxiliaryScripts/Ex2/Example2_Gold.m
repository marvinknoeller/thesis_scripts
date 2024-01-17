clear all
close all
clc
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rmpath(genpath('PrecomputedResults'));

frVec = [400 500 600 700];
for kk = 1 : 2
given_frequency = frVec(kk);
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;%.4;%.4112;
material = "gold";

radiusp = .001;    %in *100% of the wavelength
heightp = 1.5;     %in *100% of the wavelength

%-------------------------------------------------------
%% Variables for the scattering problem
if material == "silver"
    lambda_eps = geteps(lambda_point);
else
    lambda_eps = geteps_gold(lambda_point);
end
Var.eps_rel = lambda_eps;
Var.mu_rel = 1;
Var.kappa = 1;
%% Variables for the optimization

aa = 1;
bb = -aa/real(lambda_eps) * 1.;
Var.N = 5;
rho = .2;
Var.aa = aa * rho;
Var.bb = bb * rho;

Var.n = 20; %number of nodes
Var.M = 11; %number of points in each subsegment including first, last point
%% Initial guess of curve
rng(1234)
tt = linspace(0,1,Var.n);
kappa2 = 2*pi/(lambda_point*1e-6);
height = lambda_point * 1000 * heightp *1e-9 * kappa2;%4.;
radius = lambda_point* 1000 * radiusp *1e-9 * kappa2;
x0 =  2*(rand(1,length(tt)) - 1/2) * radius;
y0 =  2*(rand(1,length(tt)) - 1/2) * radius;
z0 = height* tt - height/2;
pcurve = [x0;y0;z0];
%%
tt = linspace(0,1,Var.n);
alpha = 0*tt;
CurveOnlyMetalBFGS(Var,pcurve,1,alpha,strcat('GoldCurve',num2str(kk)));
end





