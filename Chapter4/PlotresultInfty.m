% This script plots all the figures for the infty. 
close all
if metallic == 1
    Vec = [1, 11, 16, 21, 41, length(A)];
    kappa2 = 2*pi / (.4 * 1e-6); %alpha = kappa2 in m^{-1}
    scaling = 1/kappa2 * 1e9; % * 1e9 gives the units in nm instead of m
else
    scaling = 1;
    Vec = [1 5 11 21 31 length(A)];%[1:length(A)];
end
if noise == 1
    Vec = length(A);
end
if metallic == 1
    col = '#757575';
else
    col = '#314676';
end

f=figure;
    f.Position = [680 753 334 345];
    f.Color = 'W';
for k=1:length(Vec)
    iteration = Vec(k);
    
    plot3(scaling*x0,scaling*y0,scaling*z0,'-','color',col,'LineWidth', 4)
    hold all
    plot3(scaling*ones(size(x0))*4,scaling*y0,scaling*z0,'-k', 'LineWidth', 2.5)
    plot3(scaling*x0,scaling*ones(size(y0))*4.001,scaling*z0,'-k', 'LineWidth', 2.5)
    plot3(scaling*x0,scaling*y0,scaling*ones(size(z0))*(-2),'-k', 'LineWidth', 2.5)
    X_stars = scaling*A{iteration};
    X = scaling*splinepoints(A{iteration},11);
    plot3(X(1,:),X(2,:),X(3,:),'-','LineWidth', 1,'Color','r')
    hold all
    plot3(X_stars(1,:),X_stars(2,:),X_stars(3,:), '.','MarkerSize',15, 'LineWidth', 1,'color','#A22223')
    plot3(scaling*ones(size(X(1,:)))*4,X(2,:),X(3,:),'-','color','#A22223', 'LineWidth', 2.5)
    plot3(X(1,:),scaling*ones(size(X(2,:)))*4.001,X(3,:),'-','color','#A22223', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),scaling*ones(size(X(3,:)))*(-2),'-','color','#A22223', 'LineWidth', 2.5)
    hold off
    axis equal
    axis(scaling*[-4.,4.,-4.001,4.001,-2,6])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    if metallic == 1
    text(177.6224339606142,173.4077955838911,354.1439612451049,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
    else
        text(3.566518209237215,2.622258813443522,5.935343227766963,...
        'units in m','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
    end
    ax = gca;
    if metallic == 1
    ax.XTick = [-200, -100, 0, 100, 200];
    ax.YTick = [-200, -100, 0, 100, 200];
    end
    hold off
    if k==1
        title(strcat("$\ell=$",num2str(iteration-1)," (initial guess)"),'Interpreter','Latex','FontSize',22)
    elseif k==length(Vec)
        title(strcat("$\ell=$",num2str(iteration-1)," (final result)"),'Interpreter','Latex','FontSize',22)
    else
        title(strcat("$\ell=$",num2str(iteration-1)),'Interpreter','Latex','FontSize',22)
    end
    drawnow;
    if save_results == 1
        filename = strcat('Infty/',num2str(k),name);
        print(gcf,'-djpeg',filename);
        print(gcf,'-depsc',filename);
        savefig(filename);
    end
end