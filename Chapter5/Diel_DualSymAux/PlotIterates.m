
% choice = [1,size(pp,3)];
load('DielectricCurve1Len6.mat')
choice = [1,11,31, 51, 71, size(pp,3)];
for kk = 1: length(choice)
    iteration = choice(kk);
    f=figure;
    f.Position = [680 753 334 345];
    f.Color = 'W';
    X_stars = pp(:,:,iteration);
    X = splinepoints(pp(:,:,iteration),11);
%     plot3(X(1,:),X(2,:),X(3,:),'-','LineWidth', 1,'Color','b', 'LineWidth', 2.5)
    scaling = 1;
    GenerateEllipticScattererNew(pp(:,:,iteration),Var.M,scaling*Var.aa ,scaling*Var.bb ,0*alpha,...
    RVec(:,:,iteration),SVec(:,:,iteration),TVec(:,:,iteration),'mblue',0,1)
    axis on
    camlight('right')
    hold all
    plot3(ones(size(X(1,:)))*4,X(2,:),X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),ones(size(X(2,:)))*4.001,X(3,:),'-k', 'LineWidth', 2.5)
    plot3(X(1,:),X(2,:),ones(size(X(3,:)))*(-4),'-k', 'LineWidth', 2.5)
    hold off
    axis equal
    axis([-4.,4.,-4.001,4.001,-4,4])
    set(gca,'fontsize',14)
    grid on
    set(gca,'GridAlpha', 0.2);
    set(gca,'LineWidth',2.,'TickLength',[0.025 0.04]);
    ax = gca;
    ax.XAxis.TickValues = [-4,-2,0,2,4];
    ax.ZAxis.TickValues = [-4,-2,0,2,4];
%     ax.XAxis.TickLabels(5) = {''};
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
    elseif (iteration == choice(2))
        form = '%10.0e';
        iterationnum = iteration;
        for2 = form;
        chival = num2str(round(chir(iterationnum),4),form);
    smchival = num2str(round(smooth_relax(iterationnum),4),for2);
    else
        iterationnum = iteration;
        form = 2;
        for2 = 1;
        chival = num2str(round(chir(iterationnum),2),form);
    smchival = num2str(round(smooth_relax(iterationnum),2),for2);
    end
    title(strcat("$\ell=$",num2str(iteration-1),rem),'Interpreter','Latex','FontSize',18)

    text(3.045109305040029,1.869453602902468,4.351115568611363,...
        'units in m','FontSize',16,'BackGroundColor','w',...
            'Edgecolor','k',...
            'LineStyle','-',...
            'LineWidth',1,...
            'Clipping','off')
        
    text(3.149865867740232,-6.525126819991499,-0.08944876335039, {...
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
        
        filename=strcat("ThesisPlots/1IteratesI6",num2str(kk));
        print(gcf,'-depsc',filename);
end