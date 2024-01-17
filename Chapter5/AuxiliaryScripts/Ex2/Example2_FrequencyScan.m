clc

for ii = 1:4
    fVec = [400, 500, 600, 700];
    load(strcat('SilverCurve',num2str(ii),'.mat'))
    close 
    scalVec = [1.29, 1, 0.78 , 0.6 ];
    scal = scalVec(ii);
    given_frequency = fVec(ii);
%% frequency scan
silver = 1;
if silver == 1
    lambda_points = physconst('LightSpeed')./given_frequency * 1e-6;
    lambda_eps = geteps(lambda_points);
else
    lambda_points = physconst('LightSpeed')./given_frequency * 1e-6;
    lambda_eps = geteps_gold(lambda_points);
end
if silver == 1
    name = 'silverstraightscan';
else
    name = 'goldstraightscan';
end
%% Initial guess of curve
kappa2 = 2*pi/(lambda_points*1e-6);
%%
% pcurve = pp(:,:,end);
pfinal_reality = pp(:,:,end) /kappa2;

[~,~,coefs,~,ts] = splinepoints(pcurve,Var.M);
frequency = linspace(300,801,200);
wl = physconst('LightSpeed')./frequency * 1e-6;
kfreq = 2*pi./(wl*1e-6);
Var.aa = Var.aa * scal;
Var.bb = Var.bb * scal;
%%
parfor kk = 1:length(wl)
    if silver == 1
        eps_rel = geteps(wl(kk));
    else
        eps_rel = geteps_gold(wl(kk));
    end
    pfinal_sim = pfinal_reality * kfreq(kk);
    N = ceil(max(sqrt(sum(abs(pfinal_sim).^2,1)))) + 1;
    FF = FarFieldMatrixFunction_SplineRotationRMFEpsNew(pfinal_sim,Var,0*alpha,eps_rel,N,...
        RVec(:,:,end),SVec(:,:,end),TVec(:,:,end));
    [chir_f(kk),smooth_f(kk),cint(kk),sigma_p(kk),sigma_m(kk)] = chiral(FF);
end
%%
scalNums = (1/(1e-9 * kappa2)*1e-9)^2;
scaling_of_semiaxis = kfreq / kappa2;
f=figure;
f.Position = [2230 753 334 345];
yyaxis left
plot(frequency,chir_f,'--k',frequency,smooth_f,'-b','LineWidth',2);
yyaxis right
semilogy(frequency,(1./kfreq.*sqrt(cint.*scaling_of_semiaxis.^4)).^2,':','Color',[0.5 0 0.5],'LineWidth',2)
ylim([1e-14, 1e-8])
set(gca,'Position',[0.177305383760772,0.163256032676513,0.662425160176026,0.761743967323487]);
yyaxis left
xlim([300,800]);
ylim([0,1])
ell = legend('$J_2$','$J_{\rm{HS}}$','$\left\Vert T_{D_\rho} \right\Vert_{\rm{HS}}^2$','Interpreter','Latex',...
    'Fontsize',17,'Location','NorthEast');
set(gca,'FontSize', 16);
xlabel('frequency in THz','Interpreter','Latex','Fontsize',16);
set(gca,'GridAlpha', 0.5);
grid on
ax = gca;
set(gca,'LineWidth',2,'TickLength',[0.025,0.04])
ax.XTick = [300 400 500 600 700 800];
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontSize = 16;
set(gca,'GridAlpha', 0.5);
ylabel('normalized chirality measures','Interpreter','Latex','Fontsize',14)
yyaxis right
if silver == 1
    ylabel('$\left\Vert T_{D_\rho} \right\Vert_{\rm{HS}}^2$','Interpreter','Latex','Fontsize',16);%,'Position',[845.4117636014432,50.000019073486328,-1])
else
    ylabel('$\left\Vert T_{D_\rho} \right\Vert_{\rm{HS}}^2$','Interpreter','Latex','Fontsize',16);
end
ax = gca;
ax.GridAlpha = 0.4;
ax.XMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
ax.YAxis(1).MinorTickValues = 0.1:0.2:1.6;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0.5 0 0.5];
print(gcf,'-depsc',strcat('Ex2Scan',num2str(ii)))

if exist(strcat('Example2Silver',num2str(ii))) ~= 2
    save(strcat('Example2Silver',num2str(ii),'.mat'))
end
%%
clearvars -except ii
end