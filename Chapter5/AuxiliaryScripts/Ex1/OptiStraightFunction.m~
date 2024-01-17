function  OptiStraightFunction(given_frequency,DegreeOfVecSphHar, material,tracking)
heightVec = [1/4, 1/2, 1, 2 ];
for kk = 1:4
    numpoints = 10;
    tt = linspace(0,1,numpoints);
    %%
    lambda_point = physconst('LightSpeed')./given_frequency * 1e-6;%.4;%.4112;
    
    if material == "silver"
        silver = 1; 
        lambda_eps = geteps(lambda_point);
    elseif material == "gold"
        silver = 0;
        lambda_eps = geteps_gold(lambda_point);
    else
        error("Material parameters not available. Use silver or gold.")
    end
    aa = 1;
    bb = -aa/real(lambda_eps) * 1.;
    %% straight
    kappa2 = 2*pi/(lambda_point*1e-6);
    heightMult = heightVec(kk);%1/4;
    height = lambda_point * 1000 * heightMult * 1e-9 * kappa2;
    fprintf(strcat('Height in nm:', num2str(lambda_point * 1000 * heightMult),'\n'));
    x0 = zeros(1,length(tt));
    y0 = zeros(1,length(tt));
    z0 = height * tt -height/2;
    pcurve = [x0;y0;z0];
    figure
    plot3(pcurve(1,:),pcurve(2,:),pcurve(3,:))
    drawnow
    Var.eps_rel = lambda_eps;
    Var.mu_rel = 1;
    Var.kappa = 1;
    rho = .2;
    Var.aa = aa * rho;
    Var.bb = bb * rho;
    Var.N = DegreeOfVecSphHar(kk);
    Var.M = 11;
    [~,Var.length] = CurvePenalty(pcurve,Var.M,numpoints);
%     Var.length
    [R,S,T,poncurve] = DoubleReflectionFrame(pcurve,Var.M);
    %% Initial rotation
    [~,~,coefs,~,ts] = splinepoints(pcurve,Var.M);
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,numpoints,Var.M);
    rng(345)    %we always use the same seed
    alpha = rand(1,numpoints)-.5;
    alpha = 0.01*alpha;
    [~,~,coefsalpha,~,tsalpha] = splinepoints(alpha,Var.M,ts);
    [alpha_in_between,~,~,ttalpha] = allpoints(coefsalpha,tsalpha,numpoints,Var.M);
    % account for the initial rotation -> get the frame
    for ell = 1 : length(p_in_between)
        Rnew(:,ell) = cos(alpha_in_between(ell)) * R(:,ell) + sin(alpha_in_between(ell)) * S(:,ell);
        Snew(:,ell) = -sin(alpha_in_between(ell)) * R(:,ell) + cos(alpha_in_between(ell)) * S(:,ell);
    end
    R = Rnew;
    S = Snew;
    fun_alpha = FarFieldMatrixFunction_SplineRotationRMF(pcurve, Var,R,S,T);
    [chir_alpha,smooth_relax_alpha] = chiral(fun_alpha);
    B_Hess = eye(numpoints);
    %%
    MaxIteration = 1000;
    alphavec = zeros(1,numpoints,MaxIteration);
    alphavec(:,:,1) = 0*length(alpha);
    
    RVec(:,:,1) = R;
    SVec(:,:,1) = S;
    TVec(:,:,1) = T;
    
    %reg. parameter
    lambda3 = 5e-4;
    for ell = 1 : MaxIteration
        pen3 = nPdashbP(pcurve,Var.M,numpoints,ts,Var.length,R,S,T);
        func = - smooth_relax_alpha + lambda3 * pen3;
        if ell == 1
            parfor kk1 = 1 : length(ts)
                h = zeros(1,length(ts));
                h(kk1) = 1;
                DerFFOp = DerFarFieldMatrixFunction_SplineRotationRMF(pcurve,Var,h,R,S,T);
                Der(kk1) = derive_measure(fun_alpha,DerFFOp);
                PenAlphaCurv(kk1,1) = derThetanPdashbP(pcurve,h,Var.M,numpoints,ts,Var.length,R,S,T)
            end
            Gradf = - Der.' + lambda3 * PenAlphaCurv;
            grad_f(:,:,ell) = Gradf;
        else
            Gradf = grad_f(:,:,ell);
        end
        p_k = linsolve(B_Hess,-Gradf);
        
        % Initialization
        alphaBFGS = 1e-4;
        Dphi_0 = Gradf.'*p_k;
        lambda = 0.9;
        jj = 0;
        phi_val = func;
        while phi_val > func + alphaBFGS * lambda^jj * Dphi_0
            jj = jj+1;
            % movement => functional must be evaluated for moved alpha
            [~,~,coefsalpha,~,tsalpha] = splinepoints(lambda^jj * p_k.',Var.M,ts);
            [p_k_alpha_in_between,~,~,ttalpha] = allpoints(coefsalpha,tsalpha,numpoints,Var.M);
            % built the frame
            for ell2 = 1 : length(p_in_between)
                %in this special case : T = Tnew
                Tnew(:,ell2) = der_p(:,ell2)/norm(der_p(:,ell2));
                phi(ell2) = p_k_alpha_in_between(ell2);
                Rt(:,ell2) = cos(phi(ell2)) * R(:,ell2) + sin(phi(ell2)) * S(:,ell2);
                St(:,ell2) = -sin(phi(ell2)) * R(:,ell2) + cos(phi(ell2)) * S(:,ell2);
                Rnew(:,ell2) = dotReal(T(:,ell2),Tnew(:,ell2)) * Rt(:,ell2) ...
                    - dotReal(St(:,ell2),Tnew(:,ell2)) / (1+ dotReal(T(:,ell2),Tnew(:,ell2))) * cross(T(:,ell2),Tnew(:,ell2)) ...
                    - dotReal(Rt(:,ell2),Tnew(:,ell2)) * T(:,ell2);
                Snew(:,ell2) = dotReal(T(:,ell2),Tnew(:,ell2)) * St(:,ell2) ...
                    + dotReal(Rt(:,ell2),Tnew(:,ell2)) / (1+ dotReal(T(:,ell2),Tnew(:,ell2))) * cross(T(:,ell2),Tnew(:,ell2)) ...
                    - dotReal(St(:,ell2),Tnew(:,ell2)) * T(:,ell2);
            end
            fun_alpha = FarFieldMatrixFunction_SplineRotationRMF(pcurve, Var,Rnew,Snew,Tnew);
            % value for moved alpha
            [~,phi_val] = chiral(fun_alpha);
            % reg. term for moved alpha
            pen3 = nPdashbP(pcurve,Var.M,numpoints,ts,Var.length,Rnew,Snew,Tnew);
            % try: is it smaller than value before?
            phi_val = -phi_val + lambda3*pen3;
        end
        lambda_fin = lambda^jj;
        X1 = lambda_fin*p_k.';
        
        alpha = X1;
        R = Rnew;
        S = Snew;
        
        RVec(:,:,ell+1) = R;
        SVec(:,:,ell+1) = S;
        TVec(:,:,ell+1) = T;
        
        if tracking == 1
            % track the evolution of the frame
            figure(100)
            plot3(p_in_between(1,:),p_in_between(2,:),p_in_between(3,:))
            hold on
            quiver3(p_in_between(1,:),p_in_between(2,:),p_in_between(3,:),Rnew(1,:),Rnew(2,:),Rnew(3,:))
            hold off
            drawnow
        end
        
        fun_alpha = FarFieldMatrixFunction_SplineRotationRMF(pcurve, Var,R,S,T);
        [chir_alpha,smooth_relax_alpha] = chiral(fun_alpha);
        parfor kk1 = 1 : length(ts)
            h = zeros(1,length(ts));
            h(kk1) = 1;
            DerFFOp = DerFarFieldMatrixFunction_SplineRotationRMF(pcurve,Var,h,R,S,T);
            Der(kk1) = derive_measure(fun_alpha,DerFFOp);
            PenAlphaCurv(kk1,1) = derThetanPdashbP(pcurve,h,Var.M,numpoints,ts,Var.length,R,S,T)
        end
        grad_f(:,:,ell+1) = -Der.'+ lambda3 * PenAlphaCurv;
        alphavec(:,:,ell+1) = alpha;
        
        s_k = alphavec(:,:,ell+1);
        s_kT = s_k;
        s_k = s_k.';
        if norm(alphavec(:,:,ell+1))<1e-4   % This is the condition for which we break
            fprintf("too little movement. Break.\n");
            break
        end
        norm_g = sqrt(sum(sum(Gradf.^2,1)));
        y_k = (grad_f(:,:,ell+1) - grad_f(:,:,ell));
        y_kT = y_k.';
        if (s_kT*y_k/(s_kT*s_k)) >= 1e-5*norm_g
            B_Hess = B_Hess + y_k*y_kT/(y_kT*s_k) ...
                - B_Hess*s_k*s_kT*B_Hess/(s_kT*B_Hess*s_k);
        else
            fprintf("No Update\n");
            B_Hess = B_Hess;
        end
        if min(eig(B_Hess))<0
            fprintf("Matrix not spd\n")
        end
        disp([ell chir_alpha smooth_relax_alpha])
    end
    for ell = 1 : length(R)-1
        ang(ell) = acos(R(:,ell+1).'*R(:,ell));
    end
    pdiff = poncurve(3,2:end) - poncurve(3,1:end-1);
    avgang(kk) = sum(ang./pdiff)/length(pdiff);
    pause(.5)
    if silver == 1
        save(strcat('Straights/','HeightNo',num2str(kk),'_',num2str(given_frequency)));
    else
        save(strcat('Straights/','HeightNo',num2str(kk),'_',num2str(given_frequency),'G'));
    end

end

end

