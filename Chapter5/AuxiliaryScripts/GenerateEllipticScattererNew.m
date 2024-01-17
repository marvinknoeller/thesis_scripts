function [] = GenerateEllipticScattererNew(pcurve,M,aa,bb,alpha,norVec,bnorVec,tanVec,material,GetGeo,alphaVal)
% clc
plotpoints = 0;
thetapoints = 30;
coarsepoints = 50;
finepoints = 40;
borders = pi/4;

A1 = linspace(0,borders-borders/(finepoints/2),finepoints/2);
A2 = linspace(borders,(pi-borders)-(pi - 2*borders)/coarsepoints,coarsepoints);
A3 = linspace(pi-borders,pi + borders - 2*borders/(finepoints),finepoints);
A4 = linspace(pi+borders,2*pi - borders - (pi-2*borders)/coarsepoints,coarsepoints);
A5 = linspace(2*pi - borders, 2*pi - borders/(finepoints/2),finepoints/2);
omega = [A1 A2 A3 A4 A5];
omegapoints = length(omega);
%% Definition of the center curve
tpoints = length(pcurve);
[~,~,coefs,~,ts] = splinepoints(pcurve,M);
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,tpoints,M);
[~,~,coefsA,~,ts] = splinepoints(alpha,M,ts);
[alpha,~,~,~] = allpoints(coefsA,ts,tpoints,M);
normderder_p = sqrt(sum(derder_p.^2,1));
x0 = aa * cos(omega) ;
y0 = bb * sin(omega) ;
for ell = 1 : (tpoints-1)*M - (tpoints-2)%(tpoints-1) * (M-1) 
PosOnCurve = p_in_between(:,ell);
Res = PosOnCurve + [norVec(:,ell) bnorVec(:,ell)] *[cos(alpha(ell)) -sin(alpha(ell)); sin(alpha(ell)) cos(alpha(ell))] * [x0;y0];

x0rep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(1,:);
y0rep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(2,:);
z0rep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(3,:);
end
%%
numpoints = tpoints*(M-1) - (M - 2);

for ell = 1 : numpoints
    Rtheta(1:2,1:2,ell) = eye(2);
end


%% Generate Points
PP = [x0rep; y0rep; z0rep];
% Plot points
if plotpoints == 1
f = figure;
f.Position = [2266 423 560 420];
plot3(x0rep,y0rep,z0rep,'b.','MarkerSize',15)
axis([-1,1,-1,1,-2,height+1])
axis equal
grid on
hold on
end
%% Generate Edges
% PP(:,1) is the first point, PP(:,2) is the second one and so on.
% General : PP(:,k) is the k-th point
%
% PP(:,1:omegapoints) corr. to the first ellipsoidal cross section
% PP(:,omegapoints + 1 : 2*omegapoints) corr. to the second ellipsoidal
% cross section and so on.
% General : PP(:,(k-1) * omegapoints+1 : k*omegapoints) corr. to the k-th
% ellipsoidal cross section.
%
% There must be edges connecting each of the numpoints ellipsoidal cross
% sections and edges connecting the k-th cross section with the (k-1)-th
% and (k+1)-th cross section


% connect in cross section

for k = 1 : numpoints
    crossSec = PP(:,(k-1) * omegapoints+1 : k*omegapoints);
    % we need the first point also at the end to close the ellipse
    
    
    % A numpoints*omegapoints x 2 field containing the edge correspondence
    % according to one ellipse
    EdgesSec1( (k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
        circshift(((k-1) * omegapoints+1 : k*omegapoints)',-1)];
    
    % for plotting only
    plotSecx0 = [crossSec(1,:),crossSec(1,1)];
    plotSecy0 = [crossSec(2,:),crossSec(2,1)];
    plotSecz0 = [crossSec(3,:),crossSec(3,1)];
    if plotpoints == 1
    plot3(plotSecx0,plotSecy0,plotSecz0,'-k','LineWidth',1.5)
    end
end

% connect the cross sections (numpoints - 1 connections)
for k = 1 : numpoints - 1
    crossSec_k = PP(:,(k-1) * omegapoints+1 : k*omegapoints);
    crossSec_kp1 = PP(:,(k) * omegapoints+1 : (k+1)*omegapoints);
    
    % two (numpoints-1)*omegapoints x 2 field containing the edge
    % correspondence connecting 2 ellipses
    EdgesCon2((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
        ((k) * omegapoints+1 : (k+1)*omegapoints)'];
    %Shifting
    EdgesCon1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = circshift(EdgesCon2((k-1)*omegapoints + 1 : k*omegapoints , 1:2),-1); %%
    EdgesConCross1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [circshift(((k) * omegapoints+1 : (k+1)*omegapoints)',-1),...
        ((k-1) * omegapoints+1 : k*omegapoints)'];
    
    % for plotting only (plot the cross edges)
    if plotpoints == 1
    plotCross = circshift(crossSec_kp1,[0 -1]);
    for kk = 1 : omegapoints
        plot3([plotCross(1,kk) crossSec_k(1,kk)],[plotCross(2,kk) crossSec_k(2,kk)] ,[plotCross(3,kk) crossSec_k(3,kk)],'-k','LineWidth',1.5)
    end
    
    % for plotting only (plot the con edges)
    
    for kk = 1 : omegapoints
        plot3([crossSec_k(1,kk) crossSec_kp1(1,kk)],[crossSec_k(2,kk) crossSec_kp1(2,kk)] ,[crossSec_k(3,kk) crossSec_kp1(3,kk)],'-k','LineWidth',1.5)
    end
    end
end

%% Generate Triangles for the mantle
% what is a triangle? A triangle has 3 vectices and 3 edges. 
% The edges in EdgesSec, EdgesCon and EdgesConCross are already sorted in
% such a way that a triangle is associated with the passing of four points
% namely [A, B, C, A]. Therefore, we generate our triangles simply by
% searching for these circles. 

Triangles = [EdgesSec1(1:(numpoints-1)*omegapoints,1) , EdgesCon1(:,1) , EdgesConCross1(:,1)];

Triangles2 = [EdgesConCross1(:,2) , EdgesSec1(omegapoints + 1: end,2),EdgesCon2(:,2)];


%% parameters for the two ellipoids
cc = .1;
% thetapoints = 6; %determines the number of height levels on the ellipsoids
phi = omega;
thetalow = linspace(pi/2 + pi/2/thetapoints , pi, thetapoints);
thetahigh = linspace(pi/2 - pi/2/thetapoints, 0, thetapoints);
%% Generate Points for the lower semi ellipsoid
% exclude the last one since this one has to be added solely
xEllLow = aa * sin(thetalow(1:end-1)).*cos(phi).';
xEllLowrep = reshape(xEllLow,(thetapoints-1) * omegapoints,1).';
xEllLowrep(1,end+1) = aa * sin(thetalow(end));
yEllLow = bb * sin(thetalow(1:end-1)).*sin(phi).';
yEllLowrep = reshape(yEllLow,(thetapoints-1) * omegapoints,1).';
yEllLowrep(1,end+1) = bb * sin(thetalow(end));
zEllLow =   cc * cos(thetalow(1:end-1));
zEllLowrep = repmat(zEllLow,[omegapoints,1]);
zEllLowrep = reshape(zEllLowrep,omegapoints*(thetapoints-1),1).';
zEllLowrep(1,end+1) =   cc * cos(thetalow(end));

for ell = 1 : thetapoints%(tpoints-1) * (M-1) 
PosOnCurve = p_in_between(:,1);
tan = tanVec(:,1);
nor = norVec(:,1);
bnor = bnorVec(:,1);

if ell < thetapoints
    Res = PosOnCurve + tan*zEllLowrep((ell-1)*omegapoints + 1) ...
        + [nor bnor] *[cos(alpha(1)) -sin(alpha(1)); sin(alpha(1)) cos(alpha(1))] * ...
        [xEllLowrep(1,(ell-1)*omegapoints + 1 : ell*omegapoints);yEllLowrep(1,(ell-1)*omegapoints + 1 : ell*omegapoints)];
        
    xEllLowrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(1,:);
    yEllLowrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(2,:);
    zEllLowrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(3,:);
else
    Res = PosOnCurve + tan*zEllLowrep((ell-1)*omegapoints + 1) + [nor bnor] *[cos(alpha(1)) -sin(alpha(1)); sin(alpha(1)) cos(alpha(1))] * ...
        [xEllLowrep(1,end);yEllLowrep(1,end)];
    
    xEllLowrep(1,end) = Res(1,:);
    yEllLowrep(1,end) = Res(2,:);
    zEllLowrep(1,end) = Res(3,:);
end
end


PPSemiLow = [xEllLowrep; yEllLowrep; zEllLowrep];
if plotpoints == 1
plot3(xEllLowrep,yEllLowrep,zEllLowrep,'b.','MarkerSize',15)
end
%% Generate Edges for the lower semi ellipsoid
%connect in cross section
for k = 1 : thetapoints-1 % the last one does not get connected!
    crossSecEllLow = PPSemiLow(:,(k-1) * omegapoints+1 : k*omegapoints);
    % we need the first point also at the end to close the ellipse
    
    
    % A (thetapoints-1)*omegapoints x 2 field containing the edge correspondence
    % according to one ellipse
    EdgesSecEllLow( (k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
        circshift(((k-1) * omegapoints+1 : k*omegapoints)',-1)];
    
    % for plotting only
    if plotpoints == 1
    plotSecEllLowx0 = [crossSecEllLow(1,:),crossSecEllLow(1,1)];
    plotSecEllLowy0 = [crossSecEllLow(2,:),crossSecEllLow(2,1)];
    plotSecEllLowz0 = [crossSecEllLow(3,:),crossSecEllLow(3,1)];
    plot3(plotSecEllLowx0,plotSecEllLowy0,plotSecEllLowz0,'-k','LineWidth',1.5)
    end
end
% There are already size(PP,2) points. So, we have to add this number to
% the points on the lower ellipsoid now.
EdgesSecEllLow = EdgesSecEllLow + size(PP,2); 
% axis([-1,1,-1,1,-1,0.1])

% connect the cross sections (thetapoints connections) The first and the
% last connection have to be considered seperately. In this one we only
% exclude the pole-connection. The mantle - to ellispoid connection is
% considered lateron 
for k = 1 : thetapoints -1
    if k == thetapoints-1   %connection to the pole
        crossSecEllLow_k = PPSemiLow(:,(k-1) * omegapoints+1 : k*omegapoints);
        crossSecEllLow_kp1 = PPSemiLow(:,end);
        % two (thetapoints-1)*omegapoints x 2 field containing the edge
        % correspondence connecting 2 ellipses
        EdgesConSemiLow2((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
            (k*omegapoints+1)*ones(1,omegapoints)'];
        %Shifting
        EdgesConSemiLow1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = circshift(EdgesConSemiLow2((k-1)*omegapoints + 1 : k*omegapoints , 1:2),-1); %%
        EdgesConSemiLowCross1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [(k*omegapoints+1)*ones(1,omegapoints)',...
            circshift(((k-1) * omegapoints+1 : (k)*omegapoints)',-0)];
        
        % for plotting only (plot the cross edges)
        if plotpoints == 1
        plotCrossEllLow = circshift(crossSecEllLow_kp1,[0 -1]);
        for kk = 1 : omegapoints
            plot3([plotCrossEllLow(1,1) crossSecEllLow_k(1,kk)],...
                [plotCrossEllLow(2,1) crossSecEllLow_k(2,kk)] ,[plotCrossEllLow(3,1) crossSecEllLow_k(3,kk)],'-k','LineWidth',1.5)
        end
        
        % for plotting only (plot the con edges)
        
        for kk = 1 : omegapoints
            plot3([crossSecEllLow_k(1,kk) crossSecEllLow_kp1(1,1)],...
                [crossSecEllLow_k(2,kk) crossSecEllLow_kp1(2,1)] ,[crossSecEllLow_k(3,kk) crossSecEllLow_kp1(3,1)],'-k','LineWidth',1.5)
        end
        end
    else
        
        crossSecEllLow_k = PPSemiLow(:,(k-1) * omegapoints+1 : k*omegapoints);
        crossSecEllLow_kp1 = PPSemiLow(:,(k) * omegapoints+1 : (k+1)*omegapoints);
        
        % two (thetapoints-1)*omegapoints x 2 field containing the edge
        % correspondence connecting 2 ellipses
        EdgesConSemiLow2((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
            ((k) * omegapoints+1 : (k+1)*omegapoints)'];
        %Shifting
        EdgesConSemiLow1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = circshift(EdgesConSemiLow2((k-1)*omegapoints + 1 : k*omegapoints , 1:2),-1); %%
        EdgesConSemiLowCross1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [circshift(((k) * omegapoints+1 : (k+1)*omegapoints)',-1),...
            ((k-1) * omegapoints+1 : k*omegapoints)'];
        % for plotting only (plot the cross edges)
        if plotpoints == 1
        plotCrossEllLow = circshift(crossSecEllLow_kp1,[0 -1]);
        for kk = 1 : omegapoints
            plot3([plotCrossEllLow(1,kk) crossSecEllLow_k(1,kk)],...
                [plotCrossEllLow(2,kk) crossSecEllLow_k(2,kk)] ,[plotCrossEllLow(3,kk) crossSecEllLow_k(3,kk)],'-k','LineWidth',1.5)
        end
        
        % for plotting only (plot the con edges)
        
        for kk = 1 : omegapoints
            plot3([crossSecEllLow_k(1,kk) crossSecEllLow_kp1(1,kk)],...
                [crossSecEllLow_k(2,kk) crossSecEllLow_kp1(2,kk)] ,[crossSecEllLow_k(3,kk) crossSecEllLow_kp1(3,kk)],'-k','LineWidth',1.5)
        end
        end
    end
    
end
% Connect mantle and ellipsoid
EdgesConSemiLow1 = EdgesConSemiLow1 + size(PP,2);
EdgesConSemiLow2 = EdgesConSemiLow2 + size(PP,2);
EdgesConSemiLowCross1 = EdgesConSemiLowCross1 + size(PP,2);

ConnectMantleEllipseLowCross = [EdgesSecEllLow(1:omegapoints,2) , (1 : omegapoints).'];
ConnectMantleEllipseLowCon = [(1 : omegapoints).',EdgesSecEllLow(1:omegapoints,1)];
ConnectMantleEllipseLowCon2 = circshift(ConnectMantleEllipseLowCon,-1);
% for plotting only (plot the cross edges)
if plotpoints == 1
for kk = 1 : omegapoints
    if kk == omegapoints
        plot3([PPSemiLow(1,1) PP(1,kk)],...
        [PPSemiLow(2,1) PP(2,kk)] ,[PPSemiLow(3,1) PP(3,kk)],'-k','LineWidth',1.5)
    else
        plot3([PPSemiLow(1,kk+1) PP(1,kk)],...
        [PPSemiLow(2,kk+1) PP(2,kk)] ,[PPSemiLow(3,kk+1) PP(3,kk)],'-k','LineWidth',1.5)
    end
end

% for plotting only (plot the con edges)

for kk = 1 : omegapoints
    plot3([PP(1,kk) PPSemiLow(1,kk)],...
        [PP(2,kk) PPSemiLow(2,kk)] ,[PP(3,kk) PPSemiLow(3,kk)],'-k','LineWidth',1.5)
end
end
%% Generate Triangles for the lower ellipsoid
TrianglesEllLow1 = [EdgesSecEllLow(:,2), EdgesConSemiLowCross1(:,2), EdgesConSemiLow1(:,2)];
TrianglesEllLow2 = [EdgesSecEllLow(omegapoints + 1: end,1), EdgesConSemiLowCross1(1:omegapoints*(thetapoints-2),1) , ...
    EdgesConSemiLow2(1:omegapoints*(thetapoints-2),1)];

TrianglesConLow = [ConnectMantleEllipseLowCon(:,1), EdgesSecEllLow(1:omegapoints,1), ConnectMantleEllipseLowCross(:,1)];
TrianglesConLow2 = [TrianglesConLow(:,3), circshift((1:omegapoints),-1).', (1:omegapoints).'];

%% Generate Points for the upper semi ellipsoid
xEllHigh = aa * sin(thetahigh(1:end-1)).*cos(phi).';
xEllHighrep = reshape(xEllHigh,(thetapoints-1) * omegapoints,1).';
xEllHighrep(1,end+1) = aa * sin(thetahigh(end));
yEllHigh = bb * sin(thetahigh(1:end-1)).*sin(phi).';
yEllHighrep = reshape(yEllHigh,(thetapoints-1) * omegapoints,1).';
yEllHighrep(1,end+1) = bb * sin(thetahigh(end));
zEllHigh =  cc * cos(thetahigh(1:end-1));
zEllHighrep = repmat(zEllHigh,[omegapoints,1]);
zEllHighrep = reshape(zEllHighrep,omegapoints*(thetapoints-1),1).';
zEllHighrep(1,end+1) =  cc * cos(thetahigh(end));

for ell = 1 : thetapoints
PosOnCurve = p_in_between(:,end);
tan = tanVec(:,end);
nor = norVec(:,end);
bnor = bnorVec(:,end);

if ell < thetapoints
    Res = PosOnCurve + tan*zEllHighrep((ell-1)*omegapoints + 1) ...
        + [nor bnor] *[cos(alpha(end)) -sin(alpha(end)); sin(alpha(end)) cos(alpha(end))] * ... %there was a 1
        [xEllHighrep(1,(ell-1)*omegapoints + 1 : ell*omegapoints);yEllHighrep(1,(ell-1)*omegapoints + 1 : ell*omegapoints)];
    
    xEllHighrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(1,:);
    yEllHighrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(2,:);
    zEllHighrep((ell-1)*omegapoints + 1: ell*omegapoints) = Res(3,:);
else
    Res = PosOnCurve + tan*zEllHighrep((ell-1)*omegapoints + 1) + [nor bnor] *[cos(alpha(end)) -sin(alpha(end)); sin(alpha(end)) cos(alpha(end))] * ... %there was a 1
        [xEllHighrep(1,end);yEllHighrep(1,end)];
   
    
    xEllHighrep(1,end) = Res(1,:);
    yEllHighrep(1,end) = Res(2,:);
    zEllHighrep(1,end) = Res(3,:);
end
end

PPSemiHigh = [xEllHighrep; yEllHighrep; zEllHighrep];
if plotpoints == 1
plot3(xEllHighrep,yEllHighrep,zEllHighrep,'b.','MarkerSize',15)
end
%% Generate Edges for the higher semi ellipsoid
%connect in cross section
for k = 1 : thetapoints-1 % the last one does not get connected!
    crossSecEllHigh = PPSemiHigh(:,(k-1) * omegapoints+1 : k*omegapoints);
    % we need the first point also at the end to close the ellipse
    
    
    % A (thetapoints-1)*omegapoints x 2 field containing the edge correspondence
    % according to one ellipse
    EdgesSecEllHigh( (k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
        circshift(((k-1) * omegapoints+1 : k*omegapoints)',-1)];
    
    % for plotting only
    plotSecEllHighx0 = [crossSecEllHigh(1,:),crossSecEllHigh(1,1)];
    plotSecEllHighy0 = [crossSecEllHigh(2,:),crossSecEllHigh(2,1)];
    plotSecEllHighz0 = [crossSecEllHigh(3,:),crossSecEllHigh(3,1)];
    if plotpoints == 1
    plot3(plotSecEllHighx0,plotSecEllHighy0,plotSecEllHighz0,'-k','LineWidth',1.5)
    end
end
% There are already size(PP,2) points. So, we have to add this number to
% the points on the lower ellipsoid now.
EdgesSecEllHigh = EdgesSecEllHigh + size(PP,2) + size(PPSemiLow,2); 
% axis([-1,1,-1,1,5,6])

% connect the cross sections (thetapoints connections) The first and the
% last connection have to be considered seperately. In this one we only
% exclude the pole-connection. The mantle - to ellispoid connection is
% considered lateron 
for k = 1 : thetapoints -1
    if k == thetapoints-1   %connection to the pole
        crossSecEllHigh_k = PPSemiHigh(:,(k-1) * omegapoints+1 : k*omegapoints);
        crossSecEllHigh_kp1 = PPSemiHigh(:,end);
        % two (thetapoints-1)*omegapoints x 2 field containing the edge
        % correspondence connecting 2 ellipses
        EdgesConSemiHigh2((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
            (k*omegapoints+1)*ones(1,omegapoints)'];
        %Shifting
        EdgesConSemiHigh1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = circshift(EdgesConSemiHigh2((k-1)*omegapoints + 1 : k*omegapoints , 1:2),-1); %%
        EdgesConSemiHighCross1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [(k*omegapoints+1)*ones(1,omegapoints)',...
            circshift(((k-1) * omegapoints+1 : (k)*omegapoints)',-0)];
        
        % for plotting only (plot the cross edges)
        if plotpoints == 1
        plotCrossEllHigh = circshift(crossSecEllHigh_kp1,[0 -1]);
        for kk = 1 : omegapoints
            plot3([plotCrossEllHigh(1,1) crossSecEllHigh_k(1,kk)],...
                [plotCrossEllHigh(2,1) crossSecEllHigh_k(2,kk)] ,[plotCrossEllHigh(3,1) crossSecEllHigh_k(3,kk)],'-k','LineWidth',1.5)
        end
        
        % for plotting only (plot the con edges)
        
        for kk = 1 : omegapoints
            plot3([crossSecEllHigh_k(1,kk) crossSecEllHigh_kp1(1,1)],...
                [crossSecEllHigh_k(2,kk) crossSecEllHigh_kp1(2,1)] ,[crossSecEllHigh_k(3,kk) crossSecEllHigh_kp1(3,1)],'-k','LineWidth',1.5)
        end
        end
    else
        
        crossSecEllHigh_k = PPSemiHigh(:,(k-1) * omegapoints+1 : k*omegapoints);
        crossSecEllHigh_kp1 = PPSemiHigh(:,(k) * omegapoints+1 : (k+1)*omegapoints);
        
        % two (thetapoints-1)*omegapoints x 2 field containing the edge
        % correspondence connecting 2 ellipses
        EdgesConSemiHigh2((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [((k-1) * omegapoints+1 : k*omegapoints)',...
            ((k) * omegapoints+1 : (k+1)*omegapoints)'];
        %Shifting
        EdgesConSemiHigh1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = circshift(EdgesConSemiHigh2((k-1)*omegapoints + 1 : k*omegapoints , 1:2),-1); %%
        EdgesConSemiHighCross1((k-1)*omegapoints + 1 : k*omegapoints , 1:2) = [circshift(((k) * omegapoints+1 : (k+1)*omegapoints)',-1),...
            ((k-1) * omegapoints+1 : k*omegapoints)'];
        % for plotting only (plot the cross edges)
        if plotpoints == 1
        plotCrossEllHigh = circshift(crossSecEllHigh_kp1,[0 -1]);
        for kk = 1 : omegapoints
            plot3([plotCrossEllHigh(1,kk) crossSecEllHigh_k(1,kk)],...
                [plotCrossEllHigh(2,kk) crossSecEllHigh_k(2,kk)] ,[plotCrossEllHigh(3,kk) crossSecEllHigh_k(3,kk)],'-k','LineWidth',1.5)
        end
        
        % for plotting only (plot the con edges)
        
        for kk = 1 : omegapoints
            plot3([crossSecEllHigh_k(1,kk) crossSecEllHigh_kp1(1,kk)],...
                [crossSecEllHigh_k(2,kk) crossSecEllHigh_kp1(2,kk)] ,[crossSecEllHigh_k(3,kk) crossSecEllHigh_kp1(3,kk)],'-k','LineWidth',1.5)
        end
        end
    end
    
end
% Connect mantle and ellipsoid
EdgesConSemiHigh1 = EdgesConSemiHigh1 + size(PP,2) + size(PPSemiLow,2); 
EdgesConSemiHigh2 = EdgesConSemiHigh2 + size(PP,2)+ size(PPSemiLow,2); 
EdgesConSemiHighCross1 = EdgesConSemiHighCross1 + size(PP,2)+ size(PPSemiLow,2); 

ConnectMantleEllipseHighCross = [EdgesSecEllHigh(1:omegapoints,2) , (size(PP,2)-omegapoints+1 : size(PP,2)).'];
ConnectMantleEllipseHighCon = [(size(PP,2)-omegapoints+1: size(PP,2)).',EdgesSecEllHigh(1:omegapoints,1)];
ConnectMantleEllipseHighCon2 = circshift(ConnectMantleEllipseHighCon,-1);
if plotpoints == 1
% for plotting only (plot the cross edges)
for kk = 1 : omegapoints
    if kk == omegapoints
        plot3([PPSemiHigh(1,1) PP(1,end-omegapoints+kk)],...
        [PPSemiHigh(2,1) PP(2,end-omegapoints+kk)] ,[PPSemiHigh(3,1) PP(3,end-omegapoints+kk)],'-k','LineWidth',1.5)
    else
        plot3([PPSemiHigh(1,kk+1) PP(1,end-omegapoints+kk)],...
        [PPSemiHigh(2,kk+1) PP(2,end-omegapoints+kk)] ,[PPSemiHigh(3,kk+1) PP(3,end-omegapoints+kk)],'-k','LineWidth',1.5)
    end
end

% for plotting only (plot the con edges)

for kk = 1 : omegapoints
    plot3([PP(1,end-omegapoints+kk) PPSemiHigh(1,kk)],...
        [PP(2,end-omegapoints+kk) PPSemiHigh(2,kk)] ,[PP(3,end-omegapoints+kk) PPSemiHigh(3,kk)],'-k','LineWidth',1.5)
end
end
%% Generate Triangles for the Higher ellipsoid
TrianglesEllHigh1 = [EdgesSecEllHigh(:,2), EdgesConSemiHighCross1(:,2), EdgesConSemiHigh1(:,2)];
TrianglesEllHigh2 = [EdgesSecEllHigh(omegapoints + 1: end,1), EdgesConSemiHighCross1(1:omegapoints*(thetapoints-2),1) , ...
    EdgesConSemiHigh2(1:omegapoints*(thetapoints-2),1)];

TrianglesConHigh = [ConnectMantleEllipseHighCon(:,1), EdgesSecEllHigh(1:omegapoints,1), ConnectMantleEllipseHighCross(:,1)];
TrianglesConHigh2 = [TrianglesConHigh(:,3), circshift((size(PP,2) - omegapoints+1:size(PP,2)),-1).',...
    (size(PP,2) - omegapoints+1:size(PP,2)).'];

% For the upper ellipsoid we flip the association
TrianglesEllHigh1 = fliplr(TrianglesEllHigh1);
TrianglesEllHigh2 = fliplr(TrianglesEllHigh2);
TrianglesConHigh = fliplr(TrianglesConHigh);
TrianglesConHigh2 = fliplr(TrianglesConHigh2);
%% Overall Points and Triangles
PointsForMsh = [PP PPSemiLow PPSemiHigh];
TrianglesForMsh = [Triangles; Triangles2; TrianglesConLow; TrianglesConLow2; TrianglesEllLow1; TrianglesEllLow2;...
    TrianglesConHigh; TrianglesConHigh2; TrianglesEllHigh1; TrianglesEllHigh2];
TR = triangulation(TrianglesForMsh,PointsForMsh.');
if nargin >4
if material == "silver"
%     trisurf(TR,'FaceColor','#4682B4','EdgeColor','none','FaceAlpha',alphaVal)    %steelblue

    trisurf(TR,'FaceColor','#C0C0C0','EdgeColor','none','FaceAlpha',alphaVal)
elseif material == "gold"
    trisurf(TR,'FaceColor','#DAA520','EdgeColor','none','FaceAlpha',alphaVal)  %goldenrod
elseif material == "mblue"
% trisurf(TR,'FaceColor','#314676','EdgeColor','none')
trisurf(TR,'FaceColor','#4682B4','EdgeColor','none')
elseif material == "mgold"
trisurf(TR,'FaceColor','#cccc33','EdgeColor','none') %shining gold
elseif material == "lavender"
trisurf(TR,'FaceColor','#f4bbff','EdgeColor','none') % brilliant lavender
end
% axis image 
camlight('left')
% axis image 
axis off
% grid off
%% Write the .msh script
if GetGeo == 1
fileID = fopen(strcat('TestTwist','.msh'),'w'); %open file for writing
fprintf(fileID, '$MeshFormat\n');
fprintf(fileID, '2.2 0 8\n');
fprintf(fileID, '$EndMeshFormat\n');
fprintf(fileID, '$Nodes\n');
% fprintf(fileID, strcat(num2str(size(PP,2) + size(PPSemiLow,2)),'\n'));
fprintf(fileID, strcat(num2str(size(PointsForMsh,2)),'\n'));

for ll = 1 : size(PointsForMsh,2)
    fprintf(fileID, strcat(num2str(ll)," ",num2str(PointsForMsh(1,ll))," ",num2str(PointsForMsh(2,ll))," ",num2str(PointsForMsh(3,ll)),'\n'));
end
fprintf(fileID, '$EndNodes\n');
% print all the triangles
fprintf(fileID, '$Elements\n');
fprintf(fileID, strcat(num2str(size(TrianglesForMsh,1)),'\n'));
for ll = 1 : size(TrianglesForMsh,1)
    fprintf(fileID, strcat(num2str(ll)," ", "2 2 1 2", " ",num2str(TrianglesForMsh(ll,1))," ",num2str(TrianglesForMsh(ll,2))," ",num2str(TrianglesForMsh(ll,3)),'\n'));
end
fprintf(fileID, '$EndElements\n');
fclose(fileID);
end
end