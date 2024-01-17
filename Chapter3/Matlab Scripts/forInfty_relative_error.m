% Plots the convergence behavior of the perturbation formula w.r. to the
% number of subsegments
%% Preamble
clear all
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% Define the 2N(N-1) points on the unit sphere
N=10;
N2=2*N;
th=linspace(pi/N,pi/N*(N-1),N-1);
Theta=ones(N2,1)*th;
th2=linspace(0,2*pi-2*pi/N2,N2)';
Phi=th2*ones(1,N-1);
% Cartesian coordinates of the sampling points and their number
Z = [sin(Theta(:)').*cos(Phi(:)'); sin(Theta(:)').*sin(Phi(:)'); cos(Theta(:)')];
n = 30;
counter=1;
w=pi/N*sqrt(sin(Theta(:)));
rad_plot = horzcat(0.03:0.01:0.1,0.12,0.14,0.16,0.18,0.2);
radii = ["003","004","005","006","007","008","009","01","012","014","016","018","02"];%,"0125","015","0175","020","025","030","040","050"];
for number = 1:length(radii)
    name = strcat('farfields/dielectric/infty_farfields/farFieldinfty_new',radii(number),'.mat');
    load(name)
    E_Bempp= FarField;
    t=linspace(0,1,n);
    R=1;
    x0 = 2*cos(2*pi*t)./(1+sin(2*pi*t).^2);
    y0 = 4*cos(2*pi*t).*sin(2*pi*t)./(1+2*sin(2*pi*t).^2);
    z0 = 4*(t).^2;
    points=[x0;y0;z0];
    Var.roh=rad_plot(number);
    Var.omega=1.0;
    Var.eps0= 8.854187817e-12;
    Var.mu0 = 4*pi*1e-7;
    eps_rel = 1;
    mu_rel = 2.1;
    frequency = 100e6;
    Var.KA = 2*pi*frequency*sqrt(Var.eps0*Var.mu0);
    Var.eps_rel = eps_rel;
    Var.mu_rel = mu_rel;
    Var.A = [-1;1i;1+1i];
    % Direction of incidence
    Var.theta=1/sqrt(3)*[1;-1;1];
    NN = size(Z,2);
    % Compute the far field with the perturbation formula
    num_x = length(points);
    [X_in_between,ww,Pol_1,Pol_2] = SetupFarField(points,11,num_x,mu_rel,eps_rel);
    E_Asi = FarField_Pert_Maxwell_E_spline(Var,X_in_between,Pol_1,Pol_2,ww,Z);
    E_Asi=E_Asi.';
    normi(counter)=err_on_ff(w.',E_Bempp-E_Asi)/err_on_ff(w.',E_Bempp);
    counter=counter+1;
end
eoc = (log(normi(2:end))-log(normi(1:end-1)))./(log(rad_plot(2:end)) - log(rad_plot(1:end-1)));
figure
set(gcf,'Position',[689 560 431 388]);
set(gcf,'color','w')
loglog(rad_plot,normi,'-db',rad_plot,rad_plot.^2*3e-0,'--r','LineWidth',2,'MarkerSize',8)
set(gca,'FontSize',14)
set(gca, 'linewidth',1.5)
xlabel('radius $\rho$ of $D_\rho$ in $\mathrm{m}$','Interpreter','Latex','FontSize',20)
ylabel('$\mathtt{RelDiff}$','Interpreter','Latex','FontSize',20)
xlim([0.03,0.2]);
ylim([1.0000e-03 0.2495]);
xticks ([0.03 0.05 0.07 0.1 0.15 0.2 0.3 0.5])
grid on
ax = gca;
ax.GridAlpha = 0.8;
ax.MinorGridAlpha = 0.2;
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.TickLength = [0.04,0.04];
ax.XAxis.MinorTickValues = [0.04 0.06 0.08 0.09 0.125 0.175];
ax.MinorGridLineStyle = '--';
line([0.035,0.06],[0.0037,0.0037],'Color','r','LineWidth',2)
line([0.06,0.06],[0.0037,0.0108],'Color','r','LineWidth',2)
text(0.065,0.006,"slope 2",'Interpreter','Latex','FontSize',20,'BackgroundColor', 'w','EdgeColor','k'...
    ,'LineWidth',1.)
axis square
filename='Inftyplots/Plot2';
print(gcf,'-djpeg',filename);
print(gcf,'-depsc',filename);
savefig(filename);