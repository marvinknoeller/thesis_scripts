clc 
close all
clear

load('HeightNo4_400.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter
% it is : 1/(2*pi/(lambda_point*1e-6)) = 1/kappa2. This is the factor that
% we need to go back to nanoscale lambda[nm].
% In order to visualize it adequately, we multiply with 1e9.
% Example: A structure has a height of 3 in computation length. Then, e.g lambda_point = .5mu meter and
% thus in real units the structure has a height of approx 2.39e-7m.
% Thus: 1/(2*pi/(lambda_point*1e-6))*1e9 = 238.73nm.
ApproxCond = Var.aa * Var.bb * 1.29^2 <= 0.05^2

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa * 1.29,scaling*Var.bb * 1.29,0*alpha,R,S,T,'silver',0,1)

ax = gca;
fig = gcf;
fig.Position = [952 401 185 430];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
ax.CameraPosition = [-359.1801 -359.5721 2.7092e+03];
axis on
hold all
xlim([-30,30])
ylim([-30,30])
zlim([-752,752])
axis on
ax.XTick = [-30 -15 0 15 30];
ax.YTick = [-30 -15 0 15 30];

grid off
h = light;
lightangle(h,70,-40)
lightangle(h,220,-90)
lightangle(h,-150,-50) 
lightangle(h,-200,-30)
text(-7.678026392290803,72.08857475678428,671.2046224962869,...%-0.02922412731732,63.66363463573042,588.9306672398051,...
        'units in nm','FontSize',15,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

grid on
filename = "silverrod400S";
print(gcf, '-depsc', filename);


%%
clc 
close all
clear

load('HeightNo4_500.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa,scaling*Var.bb,0*alpha,R,S,T,'silver',0,1)

ApproxCond = Var.aa * Var.bb <=0.05^2
ax = gca;
fig = gcf;
fig.Position = [952 401 185 430];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
ax.CameraPosition = [-359.1801 -359.5721 2.7092e+03];
axis on
hold all
xlim([-30,30])
ylim([-30,30])
zlim([-752,752])
axis on
ax.XTick = [-30 -15 0 15 30];
ax.YTick = [-30 -15 0 15 30];
grid off
h = light;
lightangle(h,70,-40)
lightangle(h,220,-90)
lightangle(h,-100,-50)
lightangle(h,-200,-30)
text(-7.678026392290803,72.08857475678428,671.2046224962869,...%-0.02922412731732,63.66363463573042,588.9306672398051,...
        'units in nm','FontSize',15,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

grid on
filename = "silverrod500S";
print(gcf, '-depsc', filename);


%%
clc 
close all
clear

load('HeightNo4_600.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa*.78,scaling*Var.bb*.78,0*alpha,R,S,T,'silver',0,1)

ApproxCond = Var.aa * Var.bb * 0.78^2 <= 0.05^2
ax = gca;
fig = gcf;
fig.Position = [952 401 185 430];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
ax.CameraPosition = [-359.1801 -359.5721 2.7092e+03];
axis on
hold all
xlim([-30,30])
ylim([-30,30])
zlim([-752,752])
axis on
ax.XTick = [-30 -15 0 15 30];
ax.YTick = [-30 -15 0 15 30];
grid off
h = light;
lightangle(h,70,-40)
lightangle(h,220,-90)
lightangle(h,-100,-50)
lightangle(h,-200,-30)
text(-7.678026392290803,72.08857475678428,671.2046224962869,...%-0.02922412731732,63.66363463573042,588.9306672398051,...
        'units in nm','FontSize',15,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

grid on
filename = "silverrod600S";
print(gcf, '-depsc', filename);

for ell = 1 : length(R)-1
        ang(ell) = acos(R(:,ell+1).'*R(:,ell));
end
pdiff = poncurve(3,2:end) - poncurve(3,1:end-1);
avgang = sum(ang./pdiff)/length(pdiff)

%%
clc 
close all
clear

load('HeightNo4_700.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa*.6,scaling*Var.bb*.6,0*alpha,R,S,T,'silver',0,1)

ApproxCond = Var.aa * Var.bb * 0.6^2 <= 0.05^2
ax = gca;
fig = gcf;
fig.Position = [952 401 185 430];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
ax.CameraPosition = [-359.1801 -359.5721 2.7092e+03];
axis on
hold all
xlim([-30,30])
ylim([-30,30])
zlim([-752,752])
axis on
ax.XTick = [-30 -15 0 15 30];
ax.YTick = [-30 -15 0 15 30];
grid off
h = light;
lightangle(h,70,-40)
lightangle(h,220,-90)
lightangle(h,250,0)
lightangle(h,-200,-30)
text(-7.678026392290803,72.08857475678428,671.2046224962869,...%-0.02922412731732,63.66363463573042,588.9306672398051,...
        'units in nm','FontSize',15,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

grid on
filename = "silverrod700S";
print(gcf, '-depsc', filename);

for ell = 1 : length(R)-1
        ang(ell) = acos(R(:,ell+1).'*R(:,ell));
end
pdiff = poncurve(3,2:end) - poncurve(3,1:end-1);
avgang = sum(ang./pdiff)/length(pdiff)
%%
clc 
close all
clear

load('HeightNo4_400G.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa*1.1,scaling*Var.bb*1.1,0*alpha,R,S,T,'gold',0,1)

ApproxCond = Var.aa * Var.bb * 1.1^2 <= 0.05^2
camlight('left')
ax = gca;
ax.CameraPosition = [-313.7414 -315.3316 2.6689e+03];%[-287.6727 -283.7029 2.7222e+03];
fig = gcf;
fig.Position = [952 435 156 396];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 11;
axis on
hold all
zlim([-1*lambda_point*1000-2, 1*lambda_point*1000+2])

text(-21.12911861082989,50.95934396479106,786.0517542001253,...%-18.81567320308318,48.479782232871685,768.7173357234151,...
        'units in nm','FontSize',14,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
axis on
ax.XTick = [-20 -10 0 10 20];
ax.YTick = [-20 -10 0 10 20];
grid off
h = light;
lightangle(h,70,-40)
grid on
filename = "goldenrod400";
print(gcf, '-depsc', filename);


%%
clc 
close all
clear

load('HeightNo4_500G.mat')
scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

GenerateEllipticScattererNew(scaling*pcurve,Var.M,scaling*Var.aa*.75,scaling*Var.bb*.75,0*alpha,R,S,T,'gold',0,1)

ApproxCond = Var.aa * Var.bb * .75^2 <= 0.05^2
camlight('left')
ax = gca;
ax.CameraPosition = [-239.2064 -240.5967 2.1042e+03];
fig = gcf;
fig.Position = [952 435 156 396];
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 11;
axis on
hold all
zlim([-1*lambda_point*1000-2, 1*lambda_point*1000+2])

text(-15.736543662832446,43.824877299944546,610.8260860415994,...
        'units in nm','FontSize',14,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
axis on
ax.XTick = [-20 -10 0 10 20];
ax.YTick = [-20 -10 0 10 20];
grid off
h = light;
lightangle(h,70,-40)
grid on
filename = "goldenrod500";
print(gcf, '-depsc', filename);

for ell = 1 : length(R)-1
        ang(ell) = acos(R(:,ell+1).'*R(:,ell));
end
pdiff = poncurve(3,2:end) - poncurve(3,1:end-1);
avgang = sum(ang./pdiff)/length(pdiff)