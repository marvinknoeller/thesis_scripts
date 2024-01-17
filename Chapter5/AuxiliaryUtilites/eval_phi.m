function [chir,smooth_relax,val,FFM,dist,pen1,pen2,pen3] = eval_phi(Var,x,alpharoh,t,d,R,S,T)
    n = length(x);
    if nargin == 3 || nargin == 4
        curve = x;
    elseif nargin >4
        curve = x+t*d;
    else
        error("function is only defined for 2 or 4 inputs")
    end
    %if initialnor is defined

    if  exist('R','var')== 1
        FFM = FarFieldMatrixFunction_SplineRotationRMF(curve,Var,R,S,T);
    else
        FFM = FarFieldMatrixFunction_SplineRotationRMF(curve,Var,R,S,T);
    end
    [chir,smooth_relax] = chiral(FFM);
    Len = sum(sqrt(sum( (curve(:,2:end)-curve(:,1:end-1)).^2,1 )));
    dist = Len/(n-1);
    pen1 = Seg_Pen2(curve,Var.M,n,Var.length);
    pen2 = Curvature_Constraint(curve,Var.M,n,Var.length);
    [~,~,~,~,t_stuetz] = splinepoints(curve,Var.M);
    
%     pen3 = RotCurv(alpharoh,Var.M,n,t_stuetz);
    pen3 = nPdashbP(curve,Var.M,n,t_stuetz,Var.length,R,S,T);
    val = -(smooth_relax - Var.lambda1*pen1 - Var.lambda2*pen2 - Var.lambda3*pen3);
end