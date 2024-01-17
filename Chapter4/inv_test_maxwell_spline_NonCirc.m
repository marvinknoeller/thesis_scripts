function [ XX,errN,errD,errU,errG,j ] = inv_test_maxwell_spline_NonCirc(Var,X,FF,saveall,filename )
% inv_test returns all iterates of the newton-like algorithm and several
% error terms which appear in the function we aim to minimize.
%
% Input:
% Var - a structure which contains all the coefficients and parameters
% X - the initial guess for the curve
% FF - The far field we computed by means of Bem++
%
% Output:
% XX - a structure which contains all the iterates of the algorithm
% errN - the error of the far field
% errD - the error of the first constraint
% errU - the error of the second constraint
% errG - the full error
% j - it took j-1 steps to get to the final iterate
mu_rel = Var.mu_rel;
eps_rel = Var.eps_rel;
rho = Var.roh; 
N=10;
N2=2*N;
th=linspace(pi/N,pi/N*(N-1),N-1);
Theta=ones(N2,1)*th;
th2=linspace(0,2*pi-2*pi/N2,N2)';
Phi=th2*ones(1,N-1);
% Cartesian coordinates of the sampling points and their number
Z = [sin(Theta(:)').*cos(Phi(:)'); sin(Theta(:)').*sin(Phi(:)'); cos(Theta(:)')];
NN = size(Z,2);
%% store the given parameters in local variables
%number of segments
n=Var.n;
%the far field computed with Bem++
E_infty_BEM=FF;

%% further setup
% weights of different measurements with respect to the L2 boundary norm
w = pi/N*sqrt(sin(Theta(:)));
% Number of integration points in each sub line segment
M=11;
% relative noise level. e.g. 0.1 corr. to 10% noise
if Var.noise==1
    err = (1*rand(N2*(N-1),3)-.5+1i*(1*rand(N2*(N-1),3)-.5));
    err = err/err_on_ff(w.',err.')*Var.noise_level*err_on_ff(w.',E_infty_BEM.');
    E_infty_BEM=E_infty_BEM+err;
    err_on_ff(w.',err.')
end
%%
norm_E_infty_BEM = err_on_ff(w.',E_infty_BEM.');
%%
%% Choose the regularization parameters
lambda1 = Var.lambda1;
lambda2 = Var.lambda2;
%% The first iterate
XX{1} = X;
%% The iterative reconstruction algorithm
%counts the iterates
j = 1;
% check index for finding out if we are at a local minimum
counter = 1;
% these parameters corr. to the extended algorithm.
jump_over = 1;  % just for skipping the first if condition
total_movement = 1;
relative_movement = 1;
while j < 5000 %the condition counter==0 has been replaced
    if jump_over == 0 && relative_movement>5e-3  % we moved 
        %move
    elseif relative_movement<5e-3% we did not move %urspruenglich
        if errN(j-2)<errD(j-2) && errN(j-2)<errU(j-2)
            lambda1 = lambda1*0.5;
            lambda2 = lambda2*0.5;
            fprintf('Both regularization par. lowered. \n')
        elseif errN(j-2)<errD(j-2) && errN(j-2)>=errU(j-2)
            lambda1 = lambda1*0.5;
            fprintf('Curvature regularization par. lowered. \n')
        elseif errN(j-2)>=errD(j-2) && errN(j-2)<errU(j-2)
            lambda2 = lambda2*0.5;
            fprintf('Equal distribution regularization par. lowered. \n')
        else
            break
        end
    end
    jump_over = 0;
    %% Computation of the far fields, penalty terms and derivatives
    num_x = n+1;
    if j == 1
        [R,S,T] = DoubleReflectionFrame(X,M);
        RR{1} = R;
        SS{1} = S;
        TT{1} = T;
    else
        R = RXnew;
        S = SXnew;
        T = TXnew;
    end
    
%     drawnow
    [X_in_between,ww,Pol_1,Pol_2,ts] = SetupFarFieldNonCirc(X,M,num_x,mu_rel,eps_rel,R,S,T,rho,rho);
    % far field and the Jacobian
    E_infty_Pert = FarField_Pert_Maxwell_E_splineNonCirc(Var,X_in_between,Pol_1,Pol_2,ww,Z,rho,rho);
    [J,Jp1new,Jp2new] = Jacobian_FarField_maxwell_splineNonCirc(Var,X,Z,M,ts,R,S,T,rho,rho);
    % weighting corresponding to the L2-norm on the unit sphere
    gw = 1/norm_E_infty_BEM*(E_infty_Pert - E_infty_BEM).*w;
    gw = gw.';
    gw_vec = reshape(gw,3*NN,1);
    J = 1/norm_E_infty_BEM*w.*J;
    J_shape = reshape(permute(J,[2 1 3]),[3*NN,3*(n+1)]);
    % split it in real and imag part
    J1 = real(J_shape);
    J2 = imag(J_shape);
    gw1 = real(gw_vec);
    gw2 = imag(gw_vec);
    % first penalty term penalizes strong curvature
    p11 = Curv_Pen(X,M,num_x);
    % vector whose norm is the penalty term Psi_1
    p1 = p11(:);
    % Jacobian of the penalty term with respect to all points
    Jp1 = reshape(Jp1new,3*((num_x-1)*1*(M-1)+1),3*num_x);
    % second regularization term Psi_2 penalizes varying sublengths of
    % segments
    p22 = Seg_Pen_New(X,M,num_x);
    % vector whose norm is the penalty term Psi_2
    p2 = p22(:);
    % corresponding Jacobian
    Jp2 = Jp2new;
    %% Computation of the Gauss-Newton direction
    % Gauss-Newton iterate
    h = (J1'*J1 + J2'*J2  + lambda1^2*(Jp1'*Jp1) + lambda2^2*(Jp2'*Jp2))\...
        ([J1'  J2'  lambda1*Jp1' lambda2*Jp2']*[-gw1; -gw2; -lambda1*p1; -lambda2*p2]);
    % suggested movement of the end points
    H = reshape(h, 3, n+1);
    %% For plotting different penalty terms on logarithmic scale (for j=0 some terms equal zero)
    if j~=1
        errN(j-1)=sum(sum((abs(gw)).^2,1));
        errD(j-1)=lambda1^2*norm(p1)^2;
        errU(j-1)=lambda2^2*norm(p2)^2;
        errG(j-1)=sum(sum((abs(gw)).^2,1)) + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
    else
        errN(j)=sum(sum((abs(gw)).^2,1));
        errD(j)=lambda1^2*norm(p1)^2;
        errU(j)=lambda2^2*norm(p2)^2;
        errG(j)=sum(sum((abs(gw)).^2,1)) + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
    end
    %% Finding the minimum in the direction hh using golden ratio search
    golden = 2/(1 + sqrt(5));
    smax = 1;
    X1 = X;
    X4 = X + smax*H;
    X2 = X4 - golden*(X4 - X1); % X2 = X1 + (1-golden) * H;
    X3 = X1 + golden*(X4 - X1); % X3 = X1 + golden* H ;
    
    [~,~,coefsX2,~,tsX2] = splinepoints(X2,M,ts);
    [~,der_X2,~,~] = allpoints(coefsX2,tsX2,num_x,M);
    [~,~,coefsX3,~,tsX3] = splinepoints(X3,M,ts);
    [~,der_X3,~,~] = allpoints(coefsX3,tsX3,num_x,M);
    for ell = 1 : length(R)
        TX2(:,ell) = der_X2(:,ell)/norm(der_X2(:,ell));
        RX2(:,ell) = dotReal(T(:,ell),TX2(:,ell)) * R(:,ell) ...
            - dotReal(S(:,ell),TX2(:,ell)) / (1+ dotReal(T(:,ell),TX2(:,ell))) * cross(T(:,ell),TX2(:,ell)) ...
            - dotReal(R(:,ell),TX2(:,ell)) * T(:,ell);
        SX2(:,ell) = dotReal(T(:,ell),TX2(:,ell)) * S(:,ell) ...
            + dotReal(R(:,ell),TX2(:,ell)) / (1+ dotReal(T(:,ell),TX2(:,ell))) * cross(T(:,ell),TX2(:,ell)) ...
            - dotReal(S(:,ell),TX2(:,ell)) * T(:,ell);
        
        TX3(:,ell) = der_X3(:,ell)/norm(der_X3(:,ell));
        RX3(:,ell) = dotReal(T(:,ell),TX3(:,ell)) * R(:,ell) ...
            - dotReal(S(:,ell),TX3(:,ell)) / (1+ dotReal(T(:,ell),TX3(:,ell))) * cross(T(:,ell),TX3(:,ell)) ...
            - dotReal(R(:,ell),TX3(:,ell)) * T(:,ell);
        SX3(:,ell) = dotReal(T(:,ell),TX3(:,ell)) * S(:,ell) ...
            + dotReal(R(:,ell),TX3(:,ell)) / (1+ dotReal(T(:,ell),TX3(:,ell))) * cross(T(:,ell),TX3(:,ell)) ...
            - dotReal(S(:,ell),TX3(:,ell)) * T(:,ell);
        assert(abs(sqrt(sum(TX2(:,ell).^2))-1)<1e-10);
        assert(abs(sqrt(sum(TX3(:,ell).^2))-1)<1e-10);
        assert(abs(sqrt(sum(RX2(:,ell).^2))-1)<1e-10);
        assert(abs(sqrt(sum(RX3(:,ell).^2))-1)<1e-10);
        assert(abs(sqrt(sum(SX2(:,ell).^2))-1)<1e-10);
        assert(abs(sqrt(sum(SX3(:,ell).^2))-1)<1e-10);
        assert(abs(dotReal(TX2(:,ell),RX2(:,ell)))<1e-10);
        assert(abs(dotReal(TX2(:,ell),SX2(:,ell)))<1e-10);
        assert(abs(dotReal(RX2(:,ell),SX2(:,ell)))<1e-10);
        assert(abs(dotReal(TX3(:,ell),RX3(:,ell)))<1e-10);
        assert(abs(dotReal(TX3(:,ell),SX3(:,ell)))<1e-10);
        assert(abs(dotReal(RX3(:,ell),SX3(:,ell)))<1e-10);
    end
    [X2_in_between,ww2,Pol2_1,Pol2_2] = SetupFarFieldNonCirc(X2,M,num_x,mu_rel,eps_rel,RX2,SX2,TX2,rho,rho,ts);
    [X3_in_between,ww3,Pol3_1,Pol3_2] = SetupFarFieldNonCirc(X3,M,num_x,mu_rel,eps_rel,RX3,SX3,TX3,rho,rho,ts);
    E_infty_2 = FarField_Pert_Maxwell_E_splineNonCirc(Var,X2_in_between,Pol2_1,Pol2_2,ww2,Z,rho,rho);
    E_infty_3 = FarField_Pert_Maxwell_E_splineNonCirc(Var,X3_in_between,Pol3_1,Pol3_2,ww3,Z,rho,rho);
    % first penalty term
    p21 = Curv_Pen(X2,M,num_x);
    p21 = p21(:);
    % second penalty term
    p22 = Seg_Pen_New(X2,M,num_x);
    p22 = p22(:);
    % value of the functional to be minimized
    funct2 = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_2.')^2 + lambda1^2*norm(p21)^2 + lambda2^2*norm(p22)^2;
    % first penalty term
    p31 = Curv_Pen(X3,M,num_x);
    p31 = p31(:);
    % second penalty term
    p32 = Seg_Pen_New(X3,M,num_x);
    p32 = p32(:);
    % value of the functional to be minimized
    funct3 = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_3.')^2 + lambda1^2*norm(p31)^2 + lambda2^2*norm(p32)^2;
    counter = 0;
    % count how many times the farther point is chosen
    % ten divisions by the golden ratio search
    for k=1:10
        if funct2 <= funct3
            X4 = X3;
            X3 = X2;
            X2 = X4 - golden*(X4 - X1);
            Xnew = X2;
            funct3 = funct2;
            ind = 2;
        else
            X1 = X2;
            X2 = X3;
            X3 = X1 + golden*(X4 - X1);
            Xnew = X3;
            funct2 = funct3;
            ind = 3;
            counter = counter + 1;
        end
        [~,~,coefsXnew,~,tsXnew] = splinepoints(Xnew,M,ts);
        [~,der_Xnew,~,~] = allpoints(coefsXnew,tsXnew,num_x,M);
        for ell = 1 : length(R)
            TXnew(:,ell) = der_Xnew(:,ell)/norm(der_Xnew(:,ell));
            RXnew(:,ell) = dotReal(T(:,ell),TXnew(:,ell)) * R(:,ell) ...
                - dotReal(S(:,ell),TXnew(:,ell)) / (1+ dotReal(T(:,ell),TXnew(:,ell))) * cross(T(:,ell),TXnew(:,ell)) ...
                - dotReal(R(:,ell),TXnew(:,ell)) * T(:,ell);
            SXnew(:,ell) = dotReal(T(:,ell),TXnew(:,ell)) * S(:,ell) ...
                + dotReal(R(:,ell),TXnew(:,ell)) / (1+ dotReal(T(:,ell),TXnew(:,ell))) * cross(T(:,ell),TXnew(:,ell)) ...
                - dotReal(S(:,ell),TXnew(:,ell)) * T(:,ell);
%             [Rnew,Snew,Tnew] = DoubleReflectionFrame(Xnew,M);
            assert(abs(sqrt(sum(TXnew(:,ell).^2))-1)<1e-8);
            assert(abs(sqrt(sum(RXnew(:,ell).^2))-1)<1e-8);
            assert(abs(sqrt(sum(SXnew(:,ell).^2))-1)<1e-8);
            assert(abs(dotReal(TXnew(:,ell),RXnew(:,ell)))<1e-8);
            assert(abs(dotReal(TXnew(:,ell),SXnew(:,ell)))<1e-8);
            assert(abs(dotReal(RXnew(:,ell),SXnew(:,ell)))<1e-8);
        end
        [Xnew_in_between,wwnew,Polnew_1,Polnew_2] = SetupFarFieldNonCirc(Xnew,M,num_x,mu_rel,eps_rel,RXnew,SXnew,TXnew,rho,rho);
        E_infty_Pert = FarField_Pert_Maxwell_E_splineNonCirc(Var,Xnew_in_between,Polnew_1,Polnew_2,wwnew,Z,rho,rho);
        % first penalty term
        p11 = Curv_Pen(Xnew,M,num_x);
        p1 = p11(:);
        % second penalty term
        p22 = Seg_Pen_New(Xnew,M,num_x);
        p2 = p22(:);
        % value of the functional to be minimized
%         1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_Pert.')^2
        funkt = 1/norm_E_infty_BEM^2*err_on_ff(w.',E_infty_BEM.' - E_infty_Pert.')^2 + lambda1^2*norm(p1)^2 + lambda2^2*norm(p2)^2;
        if ind==2
            funct2 = funkt;
        else
            funct3 = funkt;
        end
    end
    % conservative choices for the parameter values and sum of the movement
    total_movement = sum(sqrt(sum((X - Xnew).^2)));
    relative_movement = sum(sqrt(sum((X - Xnew).^2)))/sum(sqrt(sum((X).^2)));
    X = Xnew;%Xnew; %before: X1
    % store the values
    XX{j+1} = X;
    RR{j+1} = RXnew;
    SS{j+1} = SXnew;
    TT{j+1} = TXnew;
    %% Monitoring
    % monitor the progress
    figure(41)
    X1 = splinepoints(X,11);
    plot3(X1(1,:),X1(2,:),X1(3,:), 'r-', 'LineWidth', 1)
    hold on
    plot3(X(1,:),X(2,:),X(3,:), 'r.', 'LineWidth', 1,'Markersize',16)
    hold off
    drawnow
    display([j, total_movement,relative_movement])
    j = j+1;
    
end
if saveall == 1
    save(filename)
end
end