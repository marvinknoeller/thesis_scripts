function [J_direction] = derive_farfieldNonCirc(p,h,Var,M,Z,ts,R,S,T,aa,bb)
% p corr. to spline points
% h corr. to distortion of p
% roh=Var.roh;
eps_rel = Var.eps_rel;
mu_rel = Var.mu_rel;
eps0 = Var.eps0;
theta=Var.theta;
KA=Var.KA;
A = Var.A;
num_x = length(p);
[p_in_between,w,Polp_1,Polp_2,ts,der_p] = SetupFarFieldNonCirc(p,M,num_x,mu_rel,eps_rel,R,S,T,aa,bb,ts);
[h_in_between,~,~,~,~,der_h] = SetupFarFieldNonCirc(h,M,num_x,mu_rel,eps_rel,R,S,T,aa,bb,ts);
const_mu = 2/(mu_rel+1);
const_eps = 2/(eps_rel+1);
D_mu = diag([1,const_mu,const_mu]);
D_eps = diag([1,const_eps,const_eps]);
%% Initializations
J_direction = zeros(size(Z,2),3);
tanvec = zeros(3,size(p_in_between,2));
lenvec = zeros(1,size(p_in_between,2));
der_Pol_1 = zeros(3,3,size(p_in_between,2));
der_Pol_2 = zeros(3,3,size(p_in_between,2));
%% Set up all the tensors
for ki = 1:size(p_in_between,2)
%     tan = der_p(:,ki)/(norm(der_p(:,ki)));
    tan = T(:,ki);
    tanvec(:,ki) = tan;
    lenvec(ki) = norm(der_p(:,ki));
   
    nor = R(:,ki);
    bnor = S(:,ki);
    %the matrix V, we need to diag.
    V = [tan,nor,bnor];
    V_ph = 1/(norm(der_p(:,ki)))*[dotReal(der_h(:,ki),nor)*nor+dotReal(der_h(:,ki),bnor)*bnor...
        , -dotReal(der_h(:,ki),nor)*tan, -dotReal(der_h(:,ki),bnor)*tan];
    der_Pol_1(:,:,ki) = (V_ph * D_mu*(V') + V*D_mu*(V_ph'));
    der_Pol_2(:,:,ki) = (V_ph * D_eps*(V') + V*D_eps*(V_ph'));
end

for j=1:size(Z,2)
    x_hat = Z(:,j);
    copy_x_hMat1 = repmat(cross([x_hat x_hat x_hat],eye(3)),1,1,length(p_in_between));
    copy_x_hMat2 = repmat(cross([x_hat x_hat x_hat],cross(eye(3),[x_hat x_hat x_hat])),1,1,length(p_in_between));
    % copy the measure point depending on the number of integration points
    x_hat = x_hat*ones(1,size(p_in_between,2));
    fact=(exp(1i*KA*dotReal(theta-x_hat,p_in_between)));
    copy_v1 = repmat((cross(theta,A)),1,1,length(p_in_between));
    copy_v2 = repmat(A,1,1,length(p_in_between));
    %% The first term of the derivative
    integrand1=-1i*KA*squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat1,Polp_1),copy_v1)).*(fact.*dotReal(h_in_between,x_hat));
    
    %% The second term of the integrand
    integrand2=squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat1,der_Pol_1),copy_v1)).*fact;
    
    %% The third term of the integrand
    integrand3=1i*KA*squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat1,Polp_1),copy_v1)).*(dotReal(theta,h_in_between).*fact);
    
    %% The fourth term of the integrand
    integrand4 = squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat1,Polp_1),copy_v1)).*fact.*(dotReal(tanvec,der_h))./lenvec;
    
    %% The first integral
    % the integrand
    Help1 = aa*bb*pi*(-KA^2)*(mu_rel-1)*(integrand1 + integrand2 + integrand3+ integrand4);
    
    %% The fifth term of the integrand
    integrand5 = -1i*KA*squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat2,Polp_2),copy_v2)).*(fact.*dotReal(x_hat,h_in_between));
    
    %% The sixth term of the integrand
    integrand6 = squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat2,der_Pol_2),copy_v2)).*fact;
    
    %% The seventh term of the integrand
    integrand7 = 1i*KA*squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat2,Polp_2),copy_v2)).*(fact.*dotReal(theta,h_in_between));
    
    %% The last term of the integrand
    integrand8 = squeeze(Tensor_Mat_Mult(Tensor_Mat_Mult(copy_x_hMat2,Polp_2),copy_v2)).*fact.*(dotReal(tanvec,der_h))./lenvec;
    
    %% The second integral
    Help2= -aa*bb*pi*KA^2*(1-eps_rel)*(integrand5 + integrand6 + integrand7 + integrand8);
    
    %% The full derivative
    integrand=sqrt(eps0)*(Help1+Help2);
    
    %% Save it
    J_direction(j,:) = sum([w;w;w].*integrand,2);
end



end

