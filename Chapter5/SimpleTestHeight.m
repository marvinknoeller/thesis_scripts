clear all
close all
%% technical setup -- just for programming
N = 5;
for nn=1:N
    NVec(nn^2:nn^2+2*nn)=nn;
end
% size of vectors and matrix
Q = 2*N*(N+2);
Qd2 = Q/2;

%% for one moment: define the helix
subnum = 20;
% t = linspace(0,4*pi,subnum);
t = linspace(0,1,subnum);
% c = 0.03;

kappa = 1.0;
mu_rel = 1.0;

fTHz = 700;%in THz
c = physconst('LightSpeed');
mum = c./fTHz * 1e-6;
kappa2 = 2*pi / (mum * 1e-6);
lowerh = 100; % in nm
upperh = 600; % in nm
rad = 30; %radius in nm
intrad = 4;% interior radius
height = linspace(lowerh * kappa2 * 1e-9, upperh * kappa2 * 1e-9,1000);
eps_rel = geteps(mum);%-5+.1;%-4+0.2*1i;%5;
Rad = rad * kappa2 * 1e-9;
Var.kappa = kappa;
M = 11;
Var.M = M;
Var.eps_rel = eps_rel;
Var.mu_rel = mu_rel;
% Var.roh = c;
Var.Q = Q;
Var.NVec = NVec;
Var.N = N;
% mu = 1;
Var.aa = intrad * kappa2 * 1e-9;
Var.bb = intrad * kappa2 * 1e-9;


parfor ite = 1 : length(height)
    ite/length(height) * 100;
    %     ite
    h = height(ite);
    x0 = Rad*cos(4*pi*t);
    y0 = Rad*sin(4*pi*t);
    z0 = t*h-h/2;
    p = [x0;y0;z0];
    alpha = 0*t;
    [~,~,coefs,~,t_stuetz] = splinepoints(p,M);
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,subnum,M);
    [~,~,coefsalpha,~,tsalpha] = splinepoints(alpha,M,t_stuetz);
    [alpha_in_between,~,~,ttalpha] = allpoints(coefsalpha,tsalpha,subnum,M);
    [R,S,T,~] = DoubleReflectionFrame(p,M);
    
    
%     T(:,1)
%     for ell = 1 : length(p_in_between)
%         Rnew(:,ell) = cos(alpha_in_between(ell)) * R(:,ell) + sin(alpha_in_between(ell)) * S(:,ell);
%         Snew(:,ell) = -sin(alpha_in_between(ell)) * R(:,ell) + cos(alpha_in_between(ell)) * S(:,ell);
%     end
%     R = Rnew;
%     S = Snew;
    assert(max(abs(sum(R.^2,1).^(1/2)-1))<=1e-13)
    assert(max(abs(sum(S.^2,1).^(1/2)-1))<=1e-13)
    
    h_der = derive_par_wr_height(t);
    FarFieldMatrix = FarFieldMatrixFunction_SplineRotationRMF(p,Var,R,S,T);
    FFM(:,:,ite) = FarFieldMatrix;
    [der_ffm] = derive_farfieldmatrix_SplineRotationRMFFD(p,h_der,Var,R,S,T);
    [full_der(ite)] = derive_measure(FarFieldMatrix,der_ffm);
    [chir(ite) , smooth_relax(ite),cint(ite)] = chiral(FarFieldMatrix);
%     GenerateEllipticScattererNew(p,M,Var.aa,Var.bb,alpha,R,S,T,'silver',0,1);
%     axis equal
%     drawnow
end
save('SimpleHeight')
%% 
f=figure;
f.Position = [950 300 253 345];%[2230 753 334 345];
% f.Position = [950 300 353 345];
set(gcf,'color','w');
% axes('Position',[.2 .575 .3 .32])
% box on
% ax2 = gca;
% ax2.XTick = [];
% ax2.YTick = [];
% set(gca,'linewidth',2)
% % axes('Position',[.0 .85 .2 .35])
% axes('Position',[.2 .2 .4 .65])
% box on
% plot(x2,y2)
h = 3;
alpha = 0*t;
x0 = Rad*cos(4*pi*t);
y0 = Rad*sin(4*pi*t);
z0 = t*h-h/2;
p = [x0;y0;z0];
[R,S,T,~] = DoubleReflectionFrame(p,M);
Var.aa = intrad * kappa2 * 1e-9;
Var.bb = intrad * kappa2 * 1e-9;
GenerateEllipticScattererNew(p,M,Var.aa,Var.bb,alpha,R,S,T,'silver',0,1);
axis equal
drawnow
% axis off
an = annotation("doublearrow");
an.Position = [0.754940711462451,0.223188405797101,-0.002946675818326,0.657971014492754];

an2 = annotation("doublearrow");
an2.Position = [0.494071146245059,0.176811594202899,0.153250810631697,0];
text(1.322878427227511,0.12748140208522,-0.365817406559199,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
    '$h$','FontSize',20,'BackGroundColor','w',...
    'Interpreter','Latex',...
    'Clipping','off')
% text(0.058391486936589,-0.08062425610224,-1.985063432008817,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
%     '$R=30\mathrm{nm}$','FontSize',18,'BackGroundColor','w',...
%     'Interpreter','Latex',...
%     'Clipping','off')
text(-0.263949457795906,0.174333965548346,-1.974596268208618,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
    '$R=30\mathrm{nm}$','FontSize',20,'BackGroundColor','w',...
    'Interpreter','Latex',...
    'Clipping','off')

filename = "SimpleHeightObj";
print(gcf,'-depsc',filename);

%%
f=figure;
f.Position = [2230 753 334 345];
plot(height/kappa2 * 1e9,chir,'--k',height/kappa2 * 1e9,smooth_relax...
    ,'-b','LineWidth',2);
xlim([lowerh,upperh]);
ylim([0,.6])
ell = legend('$X_{2,1}$','$X_{\rm{HS},1}$','Interpreter','Latex',...
    'Fontsize',20,'Location','NorthEast');
set(gca,'FontSize', 16);
xlabel('height $h$ in $\mathrm{nm}$','Interpreter','Latex','Fontsize',14);
set(gca,'GridAlpha', 0.5);
grid on

ax = gca;
set(gca,'LineWidth',2,'TickLength',[0.025,0.04])
% ax.XTick = 0:20:steps;
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontSize = 16;

ylabel('normalized chirality measures','Interpreter','Latex','Fontsize',14)
ax = gca;
ax.GridAlpha = 0.4;
ax.XMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
filename = "SimpleHeight";
print(gcf,'-depsc',filename);

%%
f=figure;
f.Position = [2230 753 334 345];
numapp = zeros(1,length(full_der));
xax = height;%Rad/kappa2 * 1e9;
numapp(1) = (smooth_relax(2) - smooth_relax(1)) / (xax(2)-xax(1));
for it = 2 : length(numapp)-1
    numapp(it) = (smooth_relax(it+1) - smooth_relax(it-1)) / (xax(it+1)-xax(it-1));
end
numapp(end) = (smooth_relax(end) - smooth_relax(end-1)) / (xax(end)-xax(end-1));
plot(height/kappa2 * 1e9,full_der...
    ,'-b','LineWidth',2);
xlim([lowerh,upperh]);

ell = legend("$X_{\rm{HS},1}'$",'Interpreter','Latex',...
    'Fontsize',20,'Location','NorthEast');

% ylim([0,1])
% ell = legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
%     'Fontsize',20,'Location','NorthEast');
set(gca,'FontSize', 16);
xlabel('height $h$ in $\mathrm{nm}$','Interpreter','Latex','Fontsize',14);
set(gca,'GridAlpha', 0.5);
grid on

ax = gca;
set(gca,'LineWidth',2,'TickLength',[0.025,0.04])
% ax.XTick = 0:20:steps;
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontSize = 16;

% ylabel('normalized chirality measures','Interpreter','Latex','Fontsize',14)
ax = gca;
ax.YTick = [-.5 , -.25, 0 , .25, .5];
ax.GridAlpha = 0.4;
ax.XMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";

fprintf("relative error excluding jump: \n")
numappex = numapp;
numappex(314:317) = 0;
full_derex = full_der;
full_derex(314:317) = 0;
fprintf(strcat(num2str(norm(numappex-full_derex)/norm(full_derex) * 100)," percent \n"))

filename = "SimpleHeightDer";
print(gcf,'-depsc',filename);