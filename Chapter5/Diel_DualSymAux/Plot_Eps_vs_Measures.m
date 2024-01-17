close all
clear
% clc
load('DielectricCurve1Len6.mat')
%%
% pcurve = pp(:,:,end);
pfinal_reality = pp(:,:,end);
epsvec = 1.01:0.1:50;
%%
parfor kk = 1:length(epsvec)
    eps_rel = epsvec(kk)
    pfinal_sim = pfinal_reality;
    N = ceil(max(sqrt(sum(abs(pfinal_sim).^2,1)))) + 1;
    FF = FarFieldMatrixFunction_SplineRotationRMFEpsNew(pfinal_sim,Var,0*alpha,eps_rel,N,...
        RVec(:,:,end),SVec(:,:,end),TVec(:,:,end));
    [chir_f(kk),smooth_f(kk),cint(kk),sigma_p(kk),sigma_m(kk)] = chiral(FF);
end
%%
scaling_of_semiaxis = 1;
f=figure;
f.Position = [2230 753 334 345];
plot(epsvec,chir_f,'--k',epsvec,smooth_f,'-b','LineWidth',2);
xlim([1,50]);
ylim([0,.3])
hold on
plot(5*ones(1,100),linspace(0,.5,100),'.','Color',[0.5 0 0.5],'LineWidth',2)
ell = legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
    'Fontsize',17,'Location','NorthEast');
set(gca,'FontSize', 16);
xlabel('$\varepsilon_r$','Interpreter','Latex','Fontsize',14);
set(gca,'GridAlpha', 0.5);
grid on
ax = gca;
set(gca,'LineWidth',2,'TickLength',[0.025,0.04])
ax.XTick = 0:10:50;
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontSize = 16;
set(gca,'GridAlpha', 0.5);
ylabel('normalized chirality measures','Interpreter','Latex','Fontsize',14)
ax = gca;
ax.GridAlpha = 0.4;
ax.XMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";
hold on

filename = strcat("ThesisPlots/EllvsMeas",num2str(1));
print(gcf,'-depsc',filename);