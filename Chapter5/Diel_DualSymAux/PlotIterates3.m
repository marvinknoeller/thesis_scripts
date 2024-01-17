clear

load('DielectricCurve3Len14.mat')
choice = [1,11,21,31,41, size(pp,3)];
% choice = [1,size(pp,3)];
for kk = 1: length(choice)
    iteration = choice(kk);
    f=figure;
    f.Position = [680 753 334 345];
    f.Color = 'W';
    X_stars = pp(:,:,iteration);
    X = splinepoints(pp(:,:,iteration),11);
    scaling = 1;
    GenerateEllipticScattererNew(pp(:,:,iteration),Var.M,scaling*Var.aa ,scaling*Var.bb ,0*alpha,...
    RVec(:,:,iteration),SVec(:,:,iteration),TVec(:,:,iteration),'lavender',0,1)
    axis on
    camlight('right')
    hold all
    plot3(ones(size(X(1,:)))*7,X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*7.001,X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-7),'-k', 'LineWidth', 2.5)
    hold off
    axis equal
    axis([-7.,7.,-7.001,7.001,-7,7])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    ax = gca;
    ax.XAxis.TickValues = -7:2:7;
    ax.YAxis.TickValues = -7:2:7;
    ax.ZAxis.TickValues = -7:2:7;
    hold off
    rem = "";
    if iteration == size(pp,3)
        form = 2;
        for2 = 2;
        iterationnum = iteration-1;
        chival = num2str(round(chir(iterationnum),2),form);
        smchival = num2str(round(smooth_relax(iterationnum),2),for2);
        rem = " (final result)";
    elseif (iteration == 1)
        form = '%10.0e';
        for2 = form;
        iterationnum = iteration;
        rem = " (initial curve)";
        chival = num2str(round(chir(iterationnum),4),form);
        smchival = num2str(round(smooth_relax(iterationnum),5),for2);
    else
        iterationnum = iteration;
        form = 2;
        for2 = 2;
        chival = num2str(round(chir(iterationnum),2),form);
    smchival = num2str(round(smooth_relax(iterationnum),2),for2);
    end
    title(strcat("$\ell=$",num2str(iteration-1),rem),'Interpreter','Latex','FontSize',18)
    text(5.243185371794539,2.938439325900731,8.137703327712416,...
        'units in m','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
        
    text(6.210397772149577,-10.103747126292944,-1.779889660749063, {...
        strcat('$$J_{\rm{2}}=$$',chival),...
        strcat('$$J_{\rm{HS}}=$$',smchival)},...
        'BackGroundColor','w',...
            'FontSize',16,'Interpreter','Latex',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off',...
            'HorizontalAlignment', 'right',...
            'VerticalAlignment', 'top')
        
        filename=strcat("ThesisPlots/3IteratesU14",num2str(kk));
        print(gcf,'-depsc',filename);
end