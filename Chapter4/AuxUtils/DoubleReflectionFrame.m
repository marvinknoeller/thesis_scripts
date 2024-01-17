function  [RVec, SVec, tanVec, PosOnCurve] = DoubleReflectionFrame(Pcurve, M, R0, S0 )
%   Performs the DoubleReflection method by Wang et. al. This is an
%       approximation to the rotation minimizing frame.
%   Input:
%            -Pcurve, the 3xn discretization points of the curve
%            -M, the points of each spline segment
%            -R0, initial vector R0 (if wanted)
%            -S0, initial vector S0 (if wanted)
%  Output:
%            -RVec, matrix containing all R vectors
%            -SVec, matrix containing all S vectors
%            -tanVec, matrix containing all tangent vectors
%            -PosOnCurve, the spline points
[~,~,coefs,~,ts] = splinepoints(Pcurve,M);
tpoints = size(Pcurve,2);
if tpoints == 2
    p_in_between = Pcurve;
    der_p(:,1) = p_in_between(:,2) - p_in_between(:,1);
    der_p(:,1) = der_p(:,1)/norm(der_p(:,1));
    der_p(:,2) = der_p(:,1);
else
    [p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,tpoints,M);
end
PosOnCurve = zeros(3,(tpoints-1)*M - (tpoints-2));
tanVec = zeros(3,(tpoints-1)*M - (tpoints-2));
for ell = 1 : (tpoints-1)*M - (tpoints-2)%(tpoints-1) * (M-1) 
    PosOnCurve(:,ell) = p_in_between(:,ell);
    tanVec(:,ell) = der_p(:,ell);
    tanVec(:,ell) = tanVec(:,ell)/norm(tanVec(:,ell));
end
if nargin == 2
    R0 = cross(cross(der_p(:,1),derder_p(:,1)),der_p(:,1));
    R0 = R0 / norm(R0);
    if max(isnan(R0)) == 1
        rng(123)
        R0 = rand(3,1);
        R0 = cross(R0,der_p(:,1));
        R0 = R0 / norm(R0);
%         R0 = [-1;0;0];
    end
    S0 = cross(tanVec(:,1),R0);
    S0 = S0 / norm(S0);
end
RVec = zeros(3,(tpoints-1)*M - (tpoints-2));
SVec = zeros(3,(tpoints-1)*M - (tpoints-2));

RVec(:,1) = R0;
SVec(:,1) = S0;
for ell = 1 : (tpoints-1)*M - (tpoints-2) - 1
    v1 = PosOnCurve(:,ell+1) - PosOnCurve(:,ell);
    c1 = norm(v1)^2;
    RiL = RVec(:,ell) - (2/c1) * dotReal(v1,RVec(:,ell)) * v1;
    TiL = tanVec(:,ell) - (2/c1) * dotReal(v1,tanVec(:,ell)) * v1;
    v2 = tanVec(:,ell+1) - TiL;
    c2 = norm(v2)^2;
    RVec(:,ell+1) = RiL - (2/c2)*dotReal(v2,RiL) * v2;
    SVec(:,ell+1) = cross(tanVec(:,ell+1),RVec(:,ell+1));
end

end

