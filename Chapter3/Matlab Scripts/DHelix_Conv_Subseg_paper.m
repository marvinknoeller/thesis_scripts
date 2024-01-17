% Plots the convergence behavior of the perturbation formula w.r. to the
% number of subsegments
%% Preamble
clear all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
% close all
%% Define the 2N(N-1) points on the unit sphere
N=10;
N2=2*N;
th=linspace(pi/N,pi/N*(N-1),N-1);
Theta=ones(N2,1)*th;
th2=linspace(0,2*pi-2*pi/N2,N2)';
Phi=th2*ones(1,N-1);
% Cartesian coordinates of the sampling points and their number
Z = [sin(Theta(:)').*cos(Phi(:)'); sin(Theta(:)').*sin(Phi(:)'); cos(Theta(:)')];
w=pi/N*sqrt(sin(Theta(:)));
%% Define the line
% Number of points
n = [4:1:34];
counter=1;
for number=n
    t=linspace(0,1,number);
    R=1;
    h=6;
    x0=R*cos(4*pi*t);
    y0=R*sin(4*pi*t);
    z0=t*h;
    points=[x0;y0;z0];
    
    %% The far field of the helix
    % Define the coefficients that appear in the perturbation formula
    Var.roh=0.03;
    Var.omega=1.0;
    Var.eps0= 8.854187817e-12;
    Var.mu0 = 4*pi*1e-7;
    mu_rel = 1;
    eps_rel = 2.1;
    Var.eps_rel = eps_rel;
    Var.mu_rel = mu_rel;
    frequency = 100e6;
    Var.KA = 2*pi*frequency*sqrt(Var.eps0*Var.mu0);
    Var.eps_rel = eps_rel;
    Var.mu_rel = mu_rel;
    Var.A = [-1;1i;1+1i];
    % Direction of incidence
    Var.theta=1/sqrt(3)*[1;-1;1];
    % NN is 2N(N-1), here in most cases 180
    NN = size(Z,2);
    % Compute the far field with the perturbation formula
    num_x = length(points);
    [X_in_between,ww,Pol_1,Pol_2] = SetupFarField(points,11,num_x,mu_rel,eps_rel);
    E_Asi = FarField_Pert_Maxwell_E_spline(Var,X_in_between,Pol_1,Pol_2,ww,Z);
    E_Asi=E_Asi.';
    load('farfields/dielectric/double_helix_farfields/farFieldDhelix_new003.mat')
    E_Bempp=FarField;%
    normi(counter) = err_on_ff(w.',E_Asi - E_Bempp)./err_on_ff(w.',E_Bempp);
counter=counter+1;
end
f = figure;
f.Position = [680 753 334 345];
set(gcf,'color','w')
os = 3.001;
zos = 1.001;
t=linspace(0,1,100);
R=1;
h=6;
x0=R*cos(4*pi*t);
y0=R*sin(4*pi*t);
z0=t*h;
X=[x0;y0;z0];
plot3(X(1,:),X(2,:),X(3,:),'-b','LineWidth', 4)
hold all
plot3(ones(size(X(1,:)))*4,X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),ones(size(X(2,:)))*4.001,X(3,:),'-k', 'LineWidth', 2.5)
plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-1),'-k', 'LineWidth', 2.5)
axis equal
axis([-4.,4.,-4.001,4.001,-1,7])
set(gca,'fontsize',14)
grid on
set(gca,'GridAlpha', 0.2);
set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
ax = gca;
ax.ZAxis.TickValues = [-1,1,3,5,7];
text(3.902519346254535,2.910362034708526,6.685515351846105,...
        'units in m','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
filename='Dhelixplots/Conv13';
print(gcf,'-djpeg',filename);
print(gcf,'-depsc',filename);
savefig(filename);




figure
set(gcf,'Position',[689 560 431 388]);
set(gcf,'color','w')
plot(n-1,normi,'-b','LineWidth',3)
set(gca,'FontSize',14)
xlabel('number of spline segments','FontSize',17)
ylabel('$\mathtt{RelDiff}$','Interpreter','Latex','FontSize',20)
set(gca, 'linewidth',1.5)
xticks(n(1:4:end))
yticks([0, 0.02, 0.04, 0.06, 0.08, 0.1])
xlim([4,max(n)-2])
ylim([0,0.11])


grid on
ax = gca;
ax.GridAlpha = 0.8;
ax.MinorGridAlpha = 0.2;
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.TickLength = [0.04,0.04];
ax.XAxis.MinorTickValues = [6 10 14 18 22 26 30];
ax.YAxis.MinorTickValues = [0.01 0.03 0.05 0.07 0.09];
ax.MinorGridLineStyle = '--';
ax.TickLength = [0.03,0.03];
axis square
filename='Dhelixplots/Conv23';
print(gcf,'-djpeg',filename);
print(gcf,'-depsc',filename);
savefig(filename);