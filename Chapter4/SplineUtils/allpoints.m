function [val,derval,derderval,t] = allpoints(coefs,ts,n,M)
%coefs is a tensor!
if size(coefs,3)==3
    c1 = coefs(:,:,1);
    c2 = coefs(:,:,2);
    c3 = coefs(:,:,3);
    
    val = zeros(3,(n-1)*M - (n-2));
    derval = zeros(3,(n-1)*M - (n-2));
    derderval =  zeros(3,(n-1)*M - (n-2));
    t = zeros(1,(n-1)*M - (n-2));
    for ite = 1:n-1 %loop over all segments
        xx = linspace(ts(ite),ts(ite+1),M);
        t(M*(ite-1)+1-(ite-1):M*ite-(ite)) = xx(1:end-1);
        
        [v1,dv1,ddv1] = splineeval(c1(ite,:),ts(ite),xx);
        [v2,dv2,ddv2] = splineeval(c2(ite,:),ts(ite),xx);
        [v3,dv3,ddv3] = splineeval(c3(ite,:),ts(ite),xx);
        
        val(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = v1(1:end-1);
        val(2,M*(ite-1)+1-(ite-1):M*ite-(ite)) = v2(1:end-1);
        val(3,M*(ite-1)+1-(ite-1):M*ite-(ite)) = v3(1:end-1);
        
        derval(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = dv1(1:end-1);
        derval(2,M*(ite-1)+1-(ite-1):M*ite-(ite)) = dv2(1:end-1);
        derval(3,M*(ite-1)+1-(ite-1):M*ite-(ite)) = dv3(1:end-1);
        
        derderval(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = ddv1(1:end-1);
        derderval(2,M*(ite-1)+1-(ite-1):M*ite-(ite)) = ddv2(1:end-1);
        derderval(3,M*(ite-1)+1-(ite-1):M*ite-(ite)) = ddv3(1:end-1);
    end
    t(end) = xx(end);
    
    val(1,end) = v1(end);
    val(2,end) = v2(end);
    val(3,end) = v3(end);
    
    derval(1,end) = dv1(end);
    derval(2,end) = dv2(end);
    derval(3,end) = dv3(end);
    
    derderval(1,end) = ddv1(end);
    derderval(2,end) = ddv2(end);
    derderval(3,end) = ddv3(end);
elseif size(coefs,3)==1
    c1 = coefs(:,:,1);
    
    val = zeros(1,(n-1)*M - (n-2));
    derval = zeros(1,(n-1)*M - (n-2));
    derderval =  zeros(1,(n-1)*M - (n-2));
    t = zeros(1,(n-1)*M - (n-2));
    for ite = 1:n-1 %loop over all segments
        xx = linspace(ts(ite),ts(ite+1),M);
        t(M*(ite-1)+1-(ite-1):M*ite-(ite)) = xx(1:end-1);
        
        [v1,dv1,ddv1] = splineeval(c1(ite,:),ts(ite),xx);
        
        val(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = v1(1:end-1);
        
        derval(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = dv1(1:end-1);
        
        derderval(1,M*(ite-1)+1-(ite-1):M*ite-(ite)) = ddv1(1:end-1);
    end
    t(end) = xx(end);
    
    val(1,end) = v1(end);

    
    derval(1,end) = dv1(end);

    
    derderval(1,end) = ddv1(end);

end