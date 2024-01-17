% This script sets up the whole algorithm. We choose the example we want to
% and simulate the reconstruction. The reference far field has been
% computed by means of the boundary element software Bempp. We computed
% several far fields and stored them in corresponding matlab files.
%% Preamble
clear
close all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% in case we want to use random data - reproduce results
rng(1234567)
% do you want to save a movie of the iterations?
movie=1;
% do you want to save pictures of the algorithm?
pictures=1;
% do you want to use noisy data?
noise=0;
% do you want to save the results?
save_results = 1;
% which example do you want to simulate?
% 1: torus
% 2: infinity with quadratic height
% 3: double-turn helix
Example=1;
metallic = 1;
%% Constants
% store everything in a structure
Var={};
if noise==1
    Var.noise=1;
    Var.noise_level=0.3;
else
    Var.noise=0;
end
%% INITIAL GUESS
num_x = 30;
% number of line segments
Var.n = num_x-1;
n = Var.n;
%% Regularization parameters
if metallic == 1
    Var.lambda1 = 0.09;
    Var.lambda2 = 0.6;
else
    Var.lambda1 = 0.06;
    Var.lambda2 = 0.9;
end
%% Example 1: Torus at 100Mhz, mur = 1.6, epsr = 2.5
if Example==1
    if metallic == 0
        % the coefficients
        t = linspace(0,1,num_x);
        name = 'Torus_100Mhz';
        if noise == 1
            name = strcat(name,'_noisy');
        end
        x11 = [0;2;0];%x22 = [0;-2;.1];%
        x22 = [1;2;0];
        X = x11*(1-t) + x22*t;
        Var.roh = 0.03;
        Var.eps0 = 8.854187817e-12;
        mu0 = 4*pi*1e-7;
        frequency = 100e6;
        Var.KA = 2*pi*frequency*sqrt(Var.eps0*mu0);
        Var.A = [-1;1i;1+1i];
        Var.mu_rel = 1.6;
        Var.eps_rel = 2.5;
        Var.theta=1/sqrt(3)*[1;-1;1];
        Var.roh = 0.03;
        t=linspace(0,2*pi,500);
        x0=1+(1)*cos(t);
        y0=1+(1)*sin(t);
        z0=-1+0*t;
        points=[x0;y0;z0];
        load('farfields/dielectric/ring_farfields/farFieldring_new003.mat')
        FF=FarField.';
    elseif metallic == 1
        t = linspace(0,1,num_x);
        name = 'Torus_Metal_400nm';
        if noise == 1
            name = strcat(name,'_noisy');
        end
        x11 = [0;2;0];
        x22 = [1;2;0];
        X = x11*(1-t) + x22*t;
        Var.roh = 0.03;
        Var.eps0 = 8.854187817e-12;
        mu0 = 4*pi*1e-7;
        Var.KA = 1.0;
        Var.A = [-1;1i;1+1i];
        Var.mu_rel = 1.0;
        Var.eps_rel = -4.4223 + 0.2104i;
        Var.theta=1/sqrt(3)*[1;-1;1];
        Var.roh = 0.03;
        t=linspace(0,2*pi,500);
        x0=1+(1)*cos(t);
        y0=1+(1)*sin(t);
        z0=-1+0*t;
        points=[x0;y0;z0];
        load('farfields/metallic/ring/farFieldring_new003.mat')
        FF=FarField.';
    end
end
%% Example 2: Lemniscate at 100Mhz. mur = 2.1, epsr = 1.0
if Example==2
    if metallic == 0
        % set the coefficients and load the far field
        name='Infty_100Mhz';
        if noise == 1
            name = strcat(name,'_noisy');
        end
        t = linspace(0,1,num_x);
        x11 = [2;0;0]; x22 = [2;2;0];
        X = x11*(1-t) + x22*t;
        Var.roh = 0.03;
        Var.eps0 = 8.854187817e-12;
        mu0 = 4*pi*1e-7;
        frequency = 100e6;
        %     Var.KA = 1;%
        Var.KA = 2*pi*frequency*sqrt(Var.eps0*mu0);
        Var.A = [-1;1i;1+1i];
        Var.eps_rel = 1;
        Var.mu_rel = 2.1;
        Var.theta=1/sqrt(3)*[1;-1;1];
        t = linspace(0,1,200);
        a = sqrt(2);
        x0 = a*sqrt(2)*cos(2*pi*t)./(1+sin(2*pi*t).^2);
        y0 = a*sqrt(2)*cos(2*pi*t).*sin(2*pi*t)./(1/2+sin(2*pi*t).^2);
        z0 = 4*(t).^2;
        points=[x0;y0;z0];
        load('farfields/dielectric/infty_farfields/farFieldinfty_new003.mat')
        FF = FarField.';
    elseif metallic == 1
        t = linspace(0.45,.55,num_x);
        name = 'Infty_Metal_400nm';
        if noise == 1
            name = strcat(name,'_noisy');
        end
        a = sqrt(2);
        X = [a*sqrt(2)*cos(2*pi*t)./(1+sin(2*pi*t).^2)...
            ; a*sqrt(2)*cos(2*pi*t).*sin(2*pi*t)./(1/2+sin(2*pi*t).^2); ...
            4*(t).^2;];
        Var.roh = 0.03;
        Var.eps0 = 8.854187817e-12;
        mu0 = 4*pi*1e-7;
        Var.KA = 1.0;
        Var.A = [-1;1i;1+1i];
        Var.mu_rel = 1.0;
        Var.eps_rel = -4.4223 + 0.2104i;
        Var.theta=1/sqrt(3)*[1;-1;1];
        Var.roh = 0.03;
        t=linspace(0,1,300);
        a = sqrt(2);
        x0 = a*sqrt(2)*cos(2*pi*t)./(1+sin(2*pi*t).^2);
        y0 = a*sqrt(2)*cos(2*pi*t).*sin(2*pi*t)./(1/2+sin(2*pi*t).^2);
        z0 = 4*(t).^2;
        points=[x0;y0;z0];
        load('farfields/metallic/infty/farFieldInftyMetal003.mat')
        FF=FarField.';
    end
end
%% Example 3: Double helix at 100Mhz. mur = 1.0, epsr = 2.1
if Example==3
    if metallic ==0
    % the coefficients
    name='Dhelix_100Mhz';
    if noise == 1
        name = strcat(name,'_noisy');
    end
    t = linspace(0,1,num_x);
    x11 = [0;-1;1]; x22 = [0;-2;1];
    X = x11*(1-t) + x22*t;
    Var.roh = 0.03;
    Var.eps0 = 8.854187817e-12;
    mu0 = 4*pi*1e-7;
    frequency = 100e6;
    Var.KA = 2*pi*frequency*sqrt(Var.eps0*mu0);
    Var.mu_rel = 1;
    Var.eps_rel = 2.1;

    Var.A = [-1;1i;1+1i];
    Var.theta=1/sqrt(3)*[1;-1;1];
    t = linspace(0,1,300);
    R = 1;
    h = 6;
    x0=R*cos(4*pi*t);
    y0=R*sin(4*pi*t);
    z0=h*t;
    points=[x0;y0;z0];
    load('farfields/dielectric/double_helix_farfields/farFieldDhelix_new003.mat')
    FF=FarField.';
    elseif metallic == 1
        t = linspace(0,1,num_x);
        name = 'Dhelix_Metal_400nm';
        if noise == 1
            name = strcat(name,'_noisy');
        end
        x11 = [0;0;0];%x22 = [0;-2;.1];%
        x22 = [0;.0;5];
        X = x11*(1-t) + x22*t;
        Var.roh = 0.03;
        Var.eps0 = 8.854187817e-12;
        mu0 = 4*pi*1e-7;
        Var.KA = 1.0;
        Var.A = [-1;1i;1+1i];
        Var.mu_rel = 1.0;
        Var.eps_rel = -4.4223 + 0.2104i;
        Var.theta=1/sqrt(3)*[1;-1;1];
        Var.roh = 0.03;
        t = linspace(0,1,300);
        R = 1;
        h = 6;
        x0=R*cos(4*pi*t);
        y0=R*sin(4*pi*t);
        z0=h*t;
        points=[x0;y0;z0];
        load('farfields/metallic/dhelix/farField_DHelixMetal003.mat')
        FF=FarField.';
    end
end


%% Reconstruction
A={};
A=[A,X];
arrow(1)=1;
[XX,errN,errD,errU,errG,arrnum]=inv_test_maxwell_spline_NonCirc(Var,X,FF,1,strcat(name,"all",".mat"));
A=[A,XX{2:end}];
if save_results==1
    save(strcat(name,'.mat'));
end
%% Video and plot?
if movie==1
    fig=figure(800);
    set(gcf,'color','white')
    set(fig,'DoubleBuffer','on');
    mov = VideoWriter(name,'MPEG-4');
    mov.Quality=100;
    mov.FrameRate = 10;
    open(mov);
    for k=1:length(A)
        X_stars = A{k};
        X = splinepoints(A{k},10);
        plot3(X(1,:),X(2,:),X(3,:), 'r-', 'LineWidth', 1)
        hold on
        plot3(X_stars(1,:),X_stars(2,:),X_stars(3,:), 'r*', 'LineWidth', 1)
        plot3(x0,y0,z0, 'LineWidth', 3)
        axis([min(x0)-1,max(x0)+1,min(y0)-1,max(y0)+1,min(z0)-1,max(z0)+7])
        hold off
        F = getframe(gcf);
        writeVideo(mov,F);
        drawnow;
    end
    close(mov);
end

if pictures==1
    if Example == 1
        PlotresultTorus
    elseif Example == 2
        PlotresultInfty
    elseif Example == 3
        PlotresultDHelix
    end
        
end