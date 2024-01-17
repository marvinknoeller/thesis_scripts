%% Silver, 700 THz Random radius, random height, random rotation

clear all
close all
clc
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
given_frequency = 700; %Frequency in THz
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;
material = "silver";
NumberOfSamples = 100;
% parpool(32)
for numS = 1 : NumberOfSamples
    rng('shuffle')
    radiusp = 1/3*rand(1);    %in *100% of the wavelength
    heightp = 2/3*rand(1);     %in *100% of the wavelength
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
    %this gives us the ratio.
    aa = 1;
    bb = 0.52;
    Var.N = 5;
    Var.n = 40; %number of nodes
    Var.M = 21; %number of points in each subsegment including first, last point
    %% Initial guess of curve in computational units
    rng(1234)
    tt = linspace(0,1,Var.n);
    kappa2 = 2*pi/(lambda_point*1e-6); % the true wave number corr. to the nanowire kappa2 acts as the scaling parameter!
    wavelength2 = lambda_point * 1e-6; % wave length in meter to get it in nm -> * 1e9
    heightm = wavelength2 *  heightp; % the height in meter to get it in nm -> * 1e9
    radiusm = wavelength2 *  radiusp; % the radius in meter to get it in nm -> * 1e9
    height = heightm * kappa2;
    radius = radiusm * kappa2;
    
    % the semiaxis lengths must be suff. small
    Rho = 0.05/(sqrt(aa*bb)*kappa2);
    aam = aa * Rho; % the first semiaxis length in meter to get it in nm -> * 1e9
    bbm = bb * Rho; % the second semiaxis length in meter to get it in nm -> * 1e9
    
    Var.aa = aam * kappa2; % we scale these as well
    Var.bb = bbm * kappa2; % we scale these as well
    x0 = radius*cos(8*pi*tt);
    y0 = radius*sin(8*pi*tt);
    z0 = height* tt - height/2;
    
    pcurve = [x0;y0;z0];
    
    rng('shuffle')
    alphaMult = (-3+6*rand(1));
    alphahandle =@(t)alphaMult*pi*t;
    [~,~,coefs,~,ts] = splinepoints(pcurve,Var.M);
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,Var.n,Var.M);
    alpha = alphahandle(ts);
    filename = strcat('Samples',num2str(given_frequency),'/','SampleNo',num2str(numS));
    FreeformMetalBFGS(Var,pcurve,1,alpha,filename);
    
    save(filename,'-append')
end






