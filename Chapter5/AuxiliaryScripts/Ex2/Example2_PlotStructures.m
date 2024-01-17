%% Silver 400 End
clc 
close all
clear

load('SilverCurve1.mat')
given_frequency = 400;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

% hold all
X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)

GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * 1.29,scaling*Var.bb * 1.29,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'silver',0,1)



display(['a = ', num2str(scaling*Var.aa * 1.29)]);
display(['b = ', num2str(scaling*Var.bb * 1.29)]);

ApproxCond = Var.aa * Var.bb * 1.29^2 <= 0.05^2
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off
text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
        
grid on
filename = "Ex2_1S";
print(gcf, '-depsc', filename);

%% Silver 500 End
clc 
close all
clear

load('SilverCurve2.mat')

given_frequency = 500;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)


GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * 1,scaling*Var.bb * 1,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'silver',0,1)
display(['a = ', num2str(scaling*Var.aa * 1)]);
display(['b = ', num2str(scaling*Var.bb * 1)]);
ApproxCond = Var.aa * Var.bb * 1 <= 0.05^2

% camlight('left')
ax = gca;

ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_2S";
print(gcf, '-depsc', filename);

%% Silver 600 End
clc 
close all
clear

load('SilverCurve3.mat')

given_frequency = 600;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];


X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * .78,scaling*Var.bb * .78,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'silver',0,1)
display(['a = ', num2str(scaling*Var.aa * .78)]);
display(['b = ', num2str(scaling*Var.bb * .78)]);
ApproxCond = Var.aa * Var.bb * 0.78^2 <= 0.05^2
% camlight('left')
ax = gca;

ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
        
grid on
filename = "Ex2_3S";
print(gcf, '-depsc', filename);


%% Silver 700 End
clc 
close all
clear

load('SilverCurve4.mat')

given_frequency = 700;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * .6,scaling*Var.bb * .6,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'silver',0,1)
display(['a = ', num2str(scaling*Var.aa * .6)]);
display(['b = ', num2str(scaling*Var.bb * .6)]);
ApproxCond = Var.aa * Var.bb * 0.34^2 <= 0.05^2
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
% ax.CameraPosition = 
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_4S";
print(gcf, '-depsc', filename);
%% Gold 400 End
clc 
close all
clear

load('GoldCurve1.mat')

given_frequency = 400;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * 1.1,scaling*Var.bb * 1.1,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'gold',0,1)

display(['a = ', num2str(scaling*Var.aa * 1.1)]);
display(['b = ', num2str(scaling*Var.bb * 1.1)]);
ApproxCond = Var.aa * Var.bb * 1.1^2 <= 0.05^2
% camlight('left')
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_5";
print(gcf, '-depsc', filename);

%% Gold 500 End
clc 
close all
clear

load('GoldCurve2.mat')

given_frequency = 500;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * .75,scaling*Var.bb * .75,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'gold',0,1)
display(['a = ', num2str(scaling*Var.aa * .75)]);
display(['b = ', num2str(scaling*Var.bb * .75)]);
ApproxCond = Var.aa * Var.bb * .75^2 <= 0.05^2
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off
text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_6";
print(gcf, '-depsc', filename);
%% Gold 600 End
clc 
close all
clear

load('GoldCurve3.mat')

given_frequency = 600;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * .39,scaling*Var.bb * .39,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'gold',0,1)

display(['a = ', num2str(scaling*Var.aa * .39)]);
display(['b = ', num2str(scaling*Var.bb * .39)]);
ApproxCond = Var.aa * Var.bb * .39^2 <= 0.05^2
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_7";
print(gcf, '-depsc', filename);
%% Gold 700 End
clc 
close all
clear

load('GoldCurve4.mat')

given_frequency = 700;
lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;

scaling = 1/(2*pi/(lambda_point*1e-6))*1e9;   %wavelength in mu meter

f=figure;
f.Position = [2230 753 334 345];

X_stars = pp(:,:,end)*scaling;
X = splinepoints(pp(:,:,end),Var.M)*scaling;
plot3(ones(size(X(1,:)))*(300),X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
hold all
plot3(X(1,:),ones(size(X(2,:)))*(300),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-300),'-k', 'LineWidth', 2.5)
GenerateEllipticScattererNew(scaling*pp(:,:,end),Var.M,scaling*Var.aa * .32,scaling*Var.bb * .32,0*alpha,...
    RVec(:,:,end),SVec(:,:,end),TVec(:,:,end),'gold',0,1)

display(['a = ', num2str(scaling*Var.aa * .32)]);
display(['b = ', num2str(scaling*Var.bb * .32)]);
ApproxCond = Var.aa * Var.bb * .32^2 <= 0.05^2
ax = gca;


ax.XLim = [-300, 300];
ax.YLim = [-300, 300];
ax.ZLim = [-300, 300];

ax.XTick = -300:100:300;
ax.YTick = -300:100:300;
ax.ZTick = -300:100:300;

fig = gcf;
ax.Box = 'on';
ax.LineWidth = 1.2;
ax.FontSize = 13;
axis on
hold all

axis on
grid off

text(110.9535212142619,-2.955774940652191,440.8252657874456,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')

text(333.8473403361496,-349.6567741908221,-122.9080826477475, {...
        strcat('$$J_{\rm{2}}=$$',num2str(round(chir(end),2),'%10.2f')),...
        strcat('$$J_{\rm{HS}}=$$',num2str(round(smooth_relax(end),2),'%10.2f'))},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
grid on
filename = "Ex2_8";
print(gcf, '-depsc', filename);