function [der_f] = eval_Gradf(Var,x,adder,FarFieldMatrix,R,S,T)
    n = length(x);
    curve = x;
    full_der = zeros(3,length(curve));
    [~,~,~,~,t_stuetz] = splinepoints(curve,Var.M);
    ex = exist('R','var')== 1;
    if Var.eps_rel == Var.mu_rel
        duality = 1;
    else
        duality = 0;
    end
    if  ex
        for kk = 1:length(curve)
            parfor j = 1:3
                %if initialnor
                
                [der_ffm] = derive_farfieldmatrix_SplineRotationRMFFD(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var,R,S,T); 
                [full_der(j,kk)] = derive_measure(FarFieldMatrix,der_ffm,duality);
                J1(j,kk) = Der_Seg_Pen2(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var.M,n,Var.length);
                J2(j,kk) = Der_Curvature_Constraint(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var.M,n,Var.length);
                J3(j,kk) = derPnPdashbP(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var.M,n,t_stuetz,Var.length,R,S,T)
            end
        end
    end
    der_f = -(full_der - Var.lambda1*J1 - Var.lambda2*J2 - Var.lambda3*J3);
end