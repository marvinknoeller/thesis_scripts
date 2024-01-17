% This script plots all the figures for the torus. 
close all
if metallic == 1
    Vec = [1, 11, 21, 41, 61 length(A)];
    kappa2 = 2*pi / (.4 * 1e-6); %alpha = kappa2 in m^{-1}
    scaling = 1/kappa2 * 1e9; % * 1e9 gives the units in nm instead of 
    initial_points11 = x11 * scaling
    initial_points22 = x22 * scaling
else
    scaling = 1;
    Vec = [1,6, 11, 16, 21, length(A)];
end
if noise == 1
    Vec = length(A);
end
% Vec = 1: length(A);
f=figure;
f.Position = [680 753 334 345];
f.Color = 'W';

if metallic == 1
    col = '#757575';
else
    col = '#314676';
end

for k=1:length(Vec)
    iteration = Vec(k);
    
    plot3(scaling*x0,scaling*y0,scaling*z0,'-','Color',col,'LineWidth', 4)
    hold all
    plot3(ones(size(x0))*3*scaling,y0*scaling,z0*scaling,'-k', 'LineWidth', 2.5)
    plot3(scaling*x0,ones(size(y0))*3.001*scaling,z0*scaling,'-k', 'LineWidth', 2.5)
    plot3(x0*scaling,y0*scaling,ones(size(z0))*(-3)*scaling,'-k', 'LineWidth', 2.5)
    X_stars = scaling * A{iteration};
    X = scaling * splinepoints(A{iteration},11);
    plot3(X(1,:),X(2,:),X(3,:),'-','LineWidth', 1,'Color','#A22223')
    hold all
    plot3(X_stars(1,:),X_stars(2,:),X_stars(3,:), '.','MarkerSize',15, 'LineWidth', 1,'Color','#A22223')
    plot3(ones(size(X(1,:)))*3*scaling,X(2,:),X(3,:),'-','Color','#A22223', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*3.001*scaling,X(3,:),'-','Color','#A22223', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-3*scaling),'-','Color','#A22223', 'LineWidth', 2.5)
    hold off
    axis equal
    axis(scaling*[-1.,3.,-1.001,3.001,-3,1])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    
    if metallic == 1
    text(102.3709755947921,61.03505964258693,120.3838441548869,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in nm','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
    else
      text(1.467154220946668,0.526461566439011,2.270401543639892,...%-11.321749505059643,37.301436888059925,435.9649565850114,...
        'units in m','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')  
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
        filename=strcat('Ring/',num2str(k),name);
        print(gcf,'-djpeg',filename);
        print(gcf,'-depsc',filename);
        savefig(filename);
    end
   
end