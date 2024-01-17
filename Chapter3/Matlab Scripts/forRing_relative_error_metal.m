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
alpha = 2*pi/(0.4*1e-6);
rad_plot = horzcat(0.03:0.01:0.09, 0.1,0.12,0.14,0.16,0.18,0.2);
radii = ["003","004","005","006","007","008","009","01","012","014","016","018","02"];%,"0125","015","0175","020","025","030","040","050"];
for number = 1:length(radii)
    name = strcat('farfields/metallic/ring/farFieldring_new',...
        radii(number),'.mat');
    load(name)
    E_Bempp= FarField;
    t=linspace(0,1,n);
    R=1;
    x0=R*cos(2*pi*t) +1;
    y0=R*sin(2*pi*t) +1;
    z0=0*t -1;
    points=[x0;y0;z0];
    Var.roh=rad_plot(number);
    Var.omega=1.0;
    Var.eps0= 8.854187817e-12;
    Var.mu0 = 4*pi*1e-7;
    mu_rel = 1;
    eps_rel = -4.4223 + 0.2104i;%2.5;
    frequency = 100e6;
    Var.KA = 1;%2*pi*frequency*sqrt(Var.eps0*Var.mu0);
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
rad_plot = rad_plot * 1/alpha * 1e9;
scal_to_correct_units = 1/alpha * 1e9;
figure
set(gcf,'Position',[689 560 431 388]);
set(gcf,'color','w')
loglog(rad_plot,normi,'-db',rad_plot,rad_plot.^2*6e-4,'--r','LineWidth',2,'MarkerSize',8)
hold on
% loglog(rad_plot,rad_plot*7e-3,'--m','LineWidth',2,'MarkerSize',8)
set(gca,'FontSize',14)
set(gca, 'linewidth',1.5)
xlabel('radius $\rho$ of $D_\rho$ in $\mathrm{nm}$','Interpreter','Latex','FontSize',20)
ylabel('$\mathtt{RelDiff}$','Interpreter','Latex','FontSize',20)
xlim(round(scal_to_correct_units*[0.03,0.2],2));
ylim([1.0000e-03 1]);
xticks (round(scal_to_correct_units*[0.03 0.05 0.07 0.1 0.15 0.2 0.3 0.5],2))
grid on
ax = gca;
ax.GridAlpha = 0.8;
ax.MinorGridAlpha = 0.2;
ax.XAxis.MinorTick = 'on';
ax.YAxis.MinorTick = 'on';
ax.TickLength = [0.04,0.04];
ax.XAxis.MinorTickValues = [0.04 0.06 0.08 0.09 0.125 0.175];
% ax.YAxis.MinorTickValues = logspace(-3,-1,11);%[0.005, 0.05];
ax.MinorGridLineStyle = '--';
line([2.546 4.456],[0.0038, 0.0038],'Color','r','LineWidth',2)
line([4.456, 4.456],[0.0038, 0.011],'Color','r','LineWidth',2)
text(4.7943,0.0069,0,"slope 2",'Interpreter','Latex','FontSize',20,'BackgroundColor', 'w','EdgeColor',...
    'k','LineWidth',1.)
% line([2.546, 4.456],[0.1273,0.1273],'Color','m','LineWidth',2)
% line([4.456, 4.456],[0.1273,0.2228],'Color','m','LineWidth',2)
% text(4.764692263993585,0.210939514177803,0,"slope 1",'Interpreter','Latex','FontSize',20,'BackgroundColor', 'w','EdgeColor',...
%     'k','LineWidth',1.)
axis square
filename='Ringplots/Plot1M';
print(gcf,'-djpeg',filename);
print(gcf,'-depsc',filename);
savefig(filename);