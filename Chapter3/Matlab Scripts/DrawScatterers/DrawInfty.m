%Uses pre-created mesh from gmsh
clear
close all

model = createpde;
thin = 1; % 1 for 003, 0 for 02

if thin == 1
    Infty003
    filename = 'Infty/Infty3dsthin';
else
    Infty02
    filename = 'Infty/Infty3ds';
end
nodes = msh.POS'; 

f=figure;
f.Position = [680 753 334 345];
f.Color = 'W';
TR = triangulation(msh.TRIANGLES(:,1:3),msh.POS);
trisurf(TR,'FaceColor','#4682B4','EdgeColor','none','FaceAlpha',1)

ax = gca;
axis equal
set(gca,'fontsize',14)
grid on
set(gca,'GridAlpha', 0.2);
set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
hold off
camlight left
ax.GridAlpha = 0;
axis off
print(gcf,'-depsc',filename);
savefig(filename);