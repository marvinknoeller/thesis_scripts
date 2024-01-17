function CurveOnlyDielectricBFGS(Var,pcurve,createInitialframe,alpha,filename,R,S,T)
close all
n = Var.n;
mov = 1;
for kk = 1:n
    adder(1:3,3*(kk-1)+1:3*kk,kk) = mov*eye(3);
end
B_HessP = eye(3*n);
p = pcurve;
pp(:,:,1) = p;
grad_fBFGS(:,:,1) = zeros(3,n);
lambda1 = Var.lambda1;
lambda2 = Var.lambda2;
lambda3 = Var.lambda3;
num_x = size(p,2);

if createInitialframe == 1
    [R,S,T,~] = DoubleReflectionFrame(p,Var.M);
    % Using the double reflection method to create the frame
end

RVec(:,:,1) = R;
SVec(:,:,1) = S;
TVec(:,:,1) = T;

MaxIt = 300;
for steps = 1 : MaxIt
    %% Evaluate to get the FFM
    [~,~,coefs,~,t_stuetz] = splinepoints(p,Var.M);
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);
    Pen1 = lambda1 * Seg_Pen2(p,Var.M,n,Var.length);
    Pen2 = lambda2 * Curvature_Constraint(p,Var.M,n,Var.length);
    NewVal = lambda3 * nPdashbP(p,Var.M,n,t_stuetz,Var.length,R,S,T);
    display([Pen1, Pen2, NewVal])
    % if initialnor is defined
    if steps == 1 && createInitialframe
        [~,~,coefsalpha,~,tsalpha] = splinepoints(alpha,Var.M,t_stuetz);
        [alpha_in_between,~,~,ttalpha] = allpoints(coefsalpha,tsalpha,num_x,Var.M);
        for ell = 1 : length(p_in_between)
            Rnew(:,ell) = cos(alpha_in_between(ell)) * R(:,ell) + sin(alpha_in_between(ell)) * S(:,ell);
            Snew(:,ell) = -sin(alpha_in_between(ell)) * R(:,ell) + cos(alpha_in_between(ell)) * S(:,ell);
        end
        R = Rnew;
        S = Snew;
        assert(max(abs(sum(R.^2,1).^(1/2)-1))<=1e-13)
        assert(max(abs(sum(S.^2,1).^(1/2)-1))<=1e-13)
        [chir(steps),smooth_relax(steps),func,FarFieldMatrix,~,pen1,pen2,pen3] = eval_phi(Var,p,0*alpha,0,zeros(size(p)),R,S,T);
    else
        [chir(steps),smooth_relax(steps),func,FarFieldMatrix,~,pen1,pen2,pen3] = eval_phi(Var,p,0*alpha,0,zeros(size(p)),R,S,T);
    end
    disp(strcat('Chiralitaetsmass: ' ,{' '}, num2str(chir(steps))));
    disp(strcat( 'Smooth Relaxation: ',{' '}, num2str((smooth_relax(steps)))));
    [~,~,cint] = chiral(FarFieldMatrix);
    [~,Len] = CurvePenalty(p,Var.M,num_x);
    disp(strcat( 'length of curve: ',{' '}, num2str((Len))));
    %% Store interesting values
    ValVec(steps) = smooth_relax(steps);
    Pen1Vec(steps) = lambda1*pen1;
    Pen2Vec(steps) = lambda2*pen2;
    Pen3Vec(steps) = lambda3*pen3;
    CintVec(steps) = cint;
    %%
    if steps == 1
        Gradf = eval_Gradf(Var,p,adder,FarFieldMatrix,R,S,T);
        grad_fBFGS(:,:,steps) = Gradf;
        %%
    else
        Gradf = grad_fBFGS(:,:,steps);
    end
    jac_shaped = reshape(Gradf,[],1);
    jac_shapedWR = [jac_shaped];
    p_k_shapedWR = linsolve(B_HessP,-jac_shapedWR);
    p_k_shaped = p_k_shapedWR(1:3*n);
    p_k = reshape(p_k_shaped,3,n);
    % Initialization
    alphaBFGS = 1e-4;
    phi_0 = func;
    Dphi_0 = jac_shapedWR.'*p_k_shapedWR;
    jj = 0;
    lambda = 0.9;
    if steps == 1
        pn = p + lambda^jj * p_k;
        [~,~,coefs,~,t_stuetz] = splinepoints(pn,Var.M);
        [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);
        
        
        for ell = 1 : length(p_in_between)
            Tnew(:,ell) = der_p(:,ell)/norm(der_p(:,ell));
            
            Rt(:,ell) =  R(:,ell);
            St(:,ell) =  S(:,ell);
            Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * Rt(:,ell) ...
                - dotReal(St(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                - dotReal(Rt(:,ell),Tnew(:,ell)) * T(:,ell);
            Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * St(:,ell) ...
                + dotReal(Rt(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                - dotReal(St(:,ell),Tnew(:,ell)) * T(:,ell);
        end
        [~,~,phi_val] = eval_phi(Var,p,0*alpha,lambda^jj,p_k,Rnew,Snew,Tnew);
        while phi_val > func + alphaBFGS * lambda^jj * Dphi_0
            jj = jj+1;
            pn = p + lambda^jj * p_k;
            [~,~,coefs,~,t_stuetz] = splinepoints(pn,Var.M);
            [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);

            for ell = 1 : length(p_in_between)
                Tnew(:,ell) = der_p(:,ell)/norm(der_p(:,ell));
                
                Rt(:,ell) =  R(:,ell);
                St(:,ell) =  S(:,ell);
                Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * Rt(:,ell) ...
                    - dotReal(St(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                    - dotReal(Rt(:,ell),Tnew(:,ell)) * T(:,ell);
                Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * St(:,ell) ...
                    + dotReal(Rt(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                    - dotReal(St(:,ell),Tnew(:,ell)) * T(:,ell);
            end
            [~,~,phi_val] = eval_phi(Var,p,0*alpha,lambda^jj,p_k,Rnew,Snew,Tnew);
        end
    else
        pn = p + lambda^jj * p_k;
        [~,~,coefs,~,t_stuetz] = splinepoints(pn,Var.M);
        [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);
        
        for ell = 1 : length(p_in_between)
            Tnew(:,ell) = der_p(:,ell)/norm(der_p(:,ell));
            
            Rt(:,ell) =  R(:,ell);
            St(:,ell) =  S(:,ell);
            Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * Rt(:,ell) ...
                - dotReal(St(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                - dotReal(Rt(:,ell),Tnew(:,ell)) * T(:,ell);
            Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * St(:,ell) ...
                + dotReal(Rt(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                - dotReal(St(:,ell),Tnew(:,ell)) * T(:,ell);
        end
        [~,~,phi_val] = eval_phi(Var,p,0*alpha,lambda^jj,p_k,Rnew,Snew,Tnew);
        while phi_val > func + alphaBFGS * lambda^jj * Dphi_0 && jj<200
            jj = jj+1;
            pn = p + lambda^jj * p_k;
            [~,~,coefs,~,t_stuetz] = splinepoints(pn,Var.M);
            [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);
            for ell = 1 : length(p_in_between)
                Tnew(:,ell) = der_p(:,ell)/norm(der_p(:,ell));
                
                Rt(:,ell) =  R(:,ell);
                St(:,ell) =  S(:,ell);
                Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * Rt(:,ell) ...
                    - dotReal(St(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                    - dotReal(Rt(:,ell),Tnew(:,ell)) * T(:,ell);
                Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * St(:,ell) ...
                    + dotReal(Rt(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
                    - dotReal(St(:,ell),Tnew(:,ell)) * T(:,ell);
            end
            [~,~,phi_val] = eval_phi(Var,p,0*alpha,lambda^jj,p_k,Rnew,Snew,Tnew);
        end
    end
    lambda_fin = lambda^jj;
    X1 = p + lambda_fin*p_k;    
    total_movement = sum(sqrt((sum(p - X1).^2)));
    display([steps, total_movement])
    p = X1;
    pp(:,:,steps + 1) = p;
    %% Update the frame
    [~,~,coefs,~,t_stuetz] = splinepoints(p,Var.M);
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,t_stuetz,n,Var.M);
    for ell = 1 : length(p_in_between)
        
        Tnew(:,ell) = der_p(:,ell)/norm(der_p(:,ell));
        
        Rt(:,ell) =  R(:,ell);
        St(:,ell) =  S(:,ell);
        
        Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * Rt(:,ell) ...
            - dotReal(St(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
            - dotReal(Rt(:,ell),Tnew(:,ell)) * T(:,ell);
        Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * St(:,ell) ...
            + dotReal(Rt(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
            - dotReal(St(:,ell),Tnew(:,ell)) * T(:,ell);
        assert(abs(Tnew(:,ell).'*Rnew(:,ell))<1e-12);
        assert(abs(Tnew(:,ell).'*Snew(:,ell))<1e-12);
        assert(abs(Snew(:,ell).'*Rnew(:,ell))<1e-12);
        
    end
    assert(max(abs(sum(Tnew.^2,1).^(1/2)-1))<=1e-12)
    assert(max(abs(sum(Rnew.^2,1).^(1/2)-1))<=1e-12)
    assert(max(abs(sum(Snew.^2,1).^(1/2)-1))<=1e-12)

    plottings1 = 0;
    if plottings1 == 1
        figure(39)
        subplot(1,2,1)
        plot3(p_in_between(1,:),p_in_between(2,:),p_in_between(3,:),'-*','LineWidth',2)
        view(2)
        subplot(1,2,2)
        plot3(p_in_between(1,:),p_in_between(2,:),p_in_between(3,:),'-*','LineWidth',2)
        hold on
        quiver3(p_in_between(1,:),p_in_between(2,:),p_in_between(3,:),Rnew(1,:),Rnew(2,:),Rnew(3,:),'k')
        hold off
        drawnow
    end
    R = Rnew;
    S = Snew;
    T = Tnew;
    
    RVec(:,:,steps+1) = R;
    SVec(:,:,steps+1) = S;
    TVec(:,:,steps+1) = T;
    %%
    [~,~,~,~,t_stuetz] = splinepoints(p,Var.M);
    FarFieldMatrix = FarFieldMatrixFunction_SplineRotationRMF(p,Var,R,S,T);
    
    %%
    grad_fBFGS(:,:,steps+1) = eval_Gradf(Var,p,adder,FarFieldMatrix,R,S,T);
    alphaVec(:,:,steps + 1) = alpha;
    s_k = (pp(:,:,steps+1) - pp(:,:,steps));
    s_k_shapedNoRot = reshape(s_k,[],1);
    s_k_shaped = [s_k_shapedNoRot];
    s_k_shapedT = s_k_shaped.';
    denoVec = [reshape(X1-lambda_fin*p_k,[],1)];
    if norm(s_k_shaped)/norm(denoVec)<1e-4    %1e-4
        fprintf("too little movement. Break.\n");
        break
    end
    norm_g = norm([reshape(Gradf,[],1)]);
    y_k = (grad_fBFGS(:,:,steps+1) - grad_fBFGS(:,:,steps));
    y_k_shapedNoRot = reshape(y_k,[],1);
    y_k_shaped = [y_k_shapedNoRot];
    y_k_shapedT = y_k_shaped.';
    if (s_k_shapedT*y_k_shaped/(s_k_shapedT*s_k_shaped)) >= 1e-5*norm_g
        B_HessP = B_HessP + y_k_shaped*y_k_shapedT/(y_k_shapedT*s_k_shaped) ...
            - B_HessP*s_k_shaped*s_k_shapedT*B_HessP/(s_k_shapedT*B_HessP*s_k_shaped);
    else
        fprintf("No Update\n");
        B_HessP = B_HessP;
    end
    if min(eig(B_HessP))<0
        fprintf("Matrix not spd\n")
    end
    
    plottings2 = 1;
    %%
    if plottings2 == 1
        p_spline = splinepoints(p,Var.M);
        [R,S,T,poncurve] = DoubleReflectionFrame(p,Var.M);
        figure(1)
        subplot(1,3,1)
        plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
        subplot(1,3,2)
        plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
        view(2)
        axis square
        % alpha_spline = splinepoints(alpha,Var.M);
        % subplot(1,3,3)
        % plot(tt,alpha_spline);
        figure(2)
        subplot(1,2,1)
        plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
        hold on
        quiver3(p_spline(1,:),p_spline(2,:),p_spline(3,:),R(1,:),R(2,:),R(3,:),'k')
        hold off
        subplot(1,2,2)
        plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
        % hold on
        % for kstep = 1 : length(tt)
        %     Res = [Rnew(:,kstep) Snew(:,kstep)]*[cos(thetashift) -sin(thetashift); sin(thetashift) cos(thetashift)]*...
        %         [cos(alpha_spline((kstep))) -sin(alpha_spline((kstep))); sin(alpha_spline((kstep))) cos(alpha_spline((kstep)))];
        %     RotR(:,kstep) = Res(:,1);
        % end
        % quiver3(p_spline(1,:),p_spline(2,:),p_spline(3,:),R(1,:),R(2,:),R(3,:),'k')
        hold off
        drawnow
    end
    %%
end
save(filename)
end