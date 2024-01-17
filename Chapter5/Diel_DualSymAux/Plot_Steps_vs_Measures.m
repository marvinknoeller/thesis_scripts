clear
close all
load('DielectricCurve1Len6.mat')
steps = size(pp,3)-1;
f=figure;
f.Position = [2230 753 334 345];
plot(0:steps,[chir(1:steps),chir(steps)],'--k',0:steps,[smooth_relax(1:steps),smooth_relax(steps)]...
    ,'-b','LineWidth',2);
xlim([0,steps]);
ylim([0,0.3])
ell = legend('$J_2$','$J_{\rm{HS}}$','Interpreter','Latex',...
    'Fontsize',20,'Location','NorthEast');
set(gca,'FontSize', 16);
xlabel('$\ell$','Interpreter','Latex','Fontsize',14);
set(gca,'GridAlpha', 0.5);
grid on

ax = gca;
set(gca,'LineWidth',2,'TickLength',[0.025,0.04])
ax.XTick = 0:20:steps;
ax.GridAlpha = 0.9;
ax.LineWidth = 1.2;
ax.FontSize = 16;


ylabel('normalized chirality measures','Interpreter','Latex','Fontsize',14)
ax = gca;
ax.GridAlpha = 0.4;
ax.XMinorGrid = "off";
ax.MinorGridAlpha = 0.6;
ax.MinorGridLineStyle = "--";

filename = "ThesisPlots/iter_meas1";
print(gcf,'-depsc',filename);