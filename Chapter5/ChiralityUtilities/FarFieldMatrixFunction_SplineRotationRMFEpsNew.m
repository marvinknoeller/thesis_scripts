function FarFieldMatrix = FarFieldMatrixFunction_SplineRotationRMFEpsNew(p,Var,alpharoh,eps_rel,N,R,S,T)
%% preambel
mu_rel = Var.mu_rel;
kappa = Var.kappa;
M = Var.M;
num_x = size(p,2);
for nn=1:N
    nvec(nn^2:nn^2+2*nn)=nn;
end
% size of vectors and matrix
Q = 2*N*(N+2);
Qd2 = Q/2;
% compute coefficients of diagonal matrix

aa = Var.aa;
bb = Var.bb;



const_eps1 = (aa + bb) / (aa + eps_rel*bb);%((2/((eps_rel)+1));
const_eps2 = (aa + bb) / (bb + eps_rel*aa);

const_mu1 = (aa + bb) / (aa + mu_rel*bb);%((2/((eps_rel)+1));
const_mu2 = (aa + bb) / (bb + mu_rel*aa);

[~,~,coefs,~,ts] = splinepoints(p,M);

[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
delta_t = tt(2:end) - tt(1:end-1);

[~,~,coefsalpha,~,tsalpha] = splinepoints(alpharoh,M,ts);
[alpha_in_between,~,~,ttalpha] = allpoints(coefsalpha,tsalpha,num_x,M);

n = size(p,2) - 1;
num_x = n+1;
lenV = n*M-(n-1);
for ii = 1:length(tt)
    DDreal(:,:,ii) = [real(const_eps1), 0; 0 real(const_eps2)];
    DDimag(:,:,ii) = [imag(const_eps1), 0; 0 imag(const_eps2)];
end

%%
ww = zeros(1,size(p_in_between,2));
Pol_1 = zeros(3,3,size(p_in_between,2));
Pol_2 = zeros(3,3,size(p_in_between,2));
%%
ww(1) = (delta_t(1)+delta_t(2))/6 * norm(der_p(:,1));
ww(end) = (delta_t(end)+delta_t(end-1))/6 *norm(der_p(:,end));
for ki = 2:2:(size(p_in_between,2)-1)
    ww(ki) = (delta_t(ki-1)+delta_t(ki))*4/6*norm(der_p(:,ki));   %Mid
    if ki<(size(p_in_between,2)-1)
        ww(ki+1) =  (delta_t(ki-1) + delta_t(ki) + delta_t(ki+1) + delta_t(ki+2))*1/6*norm(der_p(:,ki+1)); %End
    end
end

for ki = 1:size(p_in_between,2)

    tan = T(:,ki);
    nor = R(:,ki);
    bnor = S(:,ki);
    %%
    V = [tan,nor,bnor];
    VVV = V;
    WWW = V;
    Vmu = V;
    Pol_1(:,:,ki) = Vmu*diag([1,const_mu1,const_mu2])*(Vmu.');
    Pol_2(:,:,ki) = VVV*diag([1,DDreal(1,1,ki),DDreal(2,2,ki)])*((VVV).') + 1i*WWW*diag([0,DDimag(1,1,ki),DDimag(2,2,ki)])*(WWW.');
end
wvec = permute(repmat(ww,Q,1,Q),[1,3,2]); %Q
X1 = p_in_between(1,:);
X2 = p_in_between(2,:);
X3 = p_in_between(3,:);

lenV = (num_x-1)*M - (num_x-2);

L_inc1 = zeros(3,Qd2,lenV);
L_inc2 = zeros(3,Qd2,lenV);
R_inc1 = zeros(3,Qd2,lenV);
R_inc2 = zeros(3,Qd2,lenV);
U_scal1 = zeros(3,Qd2,lenV);
U_scal2 = zeros(3,Qd2,lenV);
B_scal1 = zeros(3,Qd2,lenV);
B_scal2 = zeros(3,Qd2,lenV);

%%
for ite=1:Qd2
    %% UPPER SCALAR PRODUCT FUNCTION
    alpha_scalU = zeros(1,Qd2);
    alpha_scalU(ite) = 1;
    beta_scalU = zeros(1,Qd2);
    scal_func1U = fields.EntireWaveField(kappa,alpha_scalU,1i*kappa*beta_scalU);
    % this is for the second term!
    scal_func2U = fields.EntireWaveField(kappa,beta_scalU,1i*kappa*alpha_scalU);
    [scal_func11U, scal_func12U, scal_func13U] = eval(scal_func1U,X1,X2,X3);
    [scal_func21U, scal_func22U, scal_func23U] = eval(scal_func2U,X1,X2,X3);
    Scal1U = zeros(3,lenV);
    Scal2U = zeros(3,lenV);
    Scal1U(1,:)=scal_func11U;
    Scal1U(2,:)=scal_func12U;
    Scal1U(3,:)=scal_func13U;
    Scal2U(1,:)=scal_func21U;
    Scal2U(2,:)=scal_func22U;
    Scal2U(3,:)=scal_func23U;
    % scaling upper left
    n_tilde = nvec(ite);
    Scal1U = 4*pi*kappa*1/((1i)^(n_tilde-1))*conj(Scal1U);
    Scal2U = (4*pi/kappa)*1/((1i)^(n_tilde-1))*conj(Scal2U);
    U_scal1(:,ite,:) = Scal1U;
    U_scal2(:,ite,:) = Scal2U;
    %% BOTTOM SCALAR PRODUCT FUNCTION
    alpha_scalB = zeros(1,Qd2);
    alpha_scalB(ite) = 1;
    beta_scalB = zeros(1,Qd2);
    scal_func1B = fields.EntireWaveField(kappa,beta_scalB,1i*kappa*alpha_scalB);    %curlMnm
    % this is for the second term!
    scal_func2B = fields.EntireWaveField(kappa,alpha_scalB,1i*kappa*beta_scalB);    %Mnm
    [scal_func11B, scal_func12B, scal_func13B] = eval(scal_func1B,X1,X2,X3);
    [scal_func21B, scal_func22B, scal_func23B] = eval(scal_func2B,X1,X2,X3);
    Scal1B = zeros(3,lenV);
    Scal2B = zeros(3,lenV);
    Scal1B(1,:)=scal_func11B;
    Scal1B(2,:)=scal_func12B;
    Scal1B(3,:)=scal_func13B;
    Scal2B(1,:)=scal_func21B;
    Scal2B(2,:)=scal_func22B;
    Scal2B(3,:)=scal_func23B;
    % scaling upper left
    n_tilde = nvec(ite);
    Scal1B = -4*pi/((1i)^(n_tilde))*conj(Scal1B);
    Scal2B = -4*pi/((1i)^(n_tilde))*conj(Scal2B);
    B_scal1(:,ite,:) = Scal1B;
    B_scal2(:,ite,:) = Scal2B;
    
    %% LEFT INCOMING FUNCTION
    alpha_scalU = zeros(1,Qd2); %could have renamed it to alpha_incL but its the same vector...
    alpha_scalU(ite) = 1;
    beta_scalU = zeros(1,Qd2);
    scal_func1U = fields.EntireWaveField(kappa,alpha_scalU,1i*kappa*beta_scalU);    %Mnm
    % this is for the second term!
    scal_func2U = fields.EntireWaveField(kappa,beta_scalU,1i*kappa*alpha_scalU);    %curlMnm
    [scal_func11U, scal_func12U, scal_func13U] = eval(scal_func1U,X1,X2,X3);
    [scal_func21U, scal_func22U, scal_func23U] = eval(scal_func2U,X1,X2,X3);
    Inc1L = zeros(3,lenV);
    Inc2L = zeros(3,lenV);
    Inc1L(1,:)=scal_func11U;
    Inc1L(2,:)=scal_func12U;
    Inc1L(3,:)=scal_func13U;
    Inc2L(1,:)=scal_func21U;
    Inc2L(2,:)=scal_func22U;
    Inc2L(3,:)=scal_func23U;
    % scaling upper left
    n_tilde = nvec(ite);
    Inc1L = -4*pi*kappa*((1i)^(n_tilde+1))*(Inc1L);
    Inc2L = (4*pi/kappa)*((1i)^(n_tilde-1))*(Inc2L);
    L_inc1(:,ite,:) = Inc1L;
    L_inc2(:,ite,:) = Inc2L;
    %% RIGHT INCOMING FUNCTION
    alpha_scalU = zeros(1,Qd2); %could have renamed it to alpha_incR but its the same vector...
    alpha_scalU(ite) = 1;
    beta_scalU = zeros(1,Qd2);
    scal_func1U = fields.EntireWaveField(kappa,beta_scalU,1i*kappa*alpha_scalU);    %curlMnm
    % this is for the second term!
    scal_func2U = fields.EntireWaveField(kappa,alpha_scalU,1i*kappa*beta_scalU);    %Mnm
    [scal_func11U, scal_func12U, scal_func13U] = eval(scal_func1U,X1,X2,X3);
    [scal_func21U, scal_func22U, scal_func23U] = eval(scal_func2U,X1,X2,X3);
    Inc1R = zeros(3,lenV);
    Inc2R = zeros(3,lenV);
    Inc1R(1,:)=scal_func11U;
    Inc1R(2,:)=scal_func12U;
    Inc1R(3,:)=scal_func13U;
    Inc2R(1,:)=scal_func21U;
    Inc2R(2,:)=scal_func22U;
    Inc2R(3,:)=scal_func23U;
    
    n_tilde = nvec(ite);
    Inc1R = -4*pi*((1i)^(n_tilde))*(Inc1R);
    Inc2R = -4*pi*((1i)^(n_tilde))*(Inc2R);
    R_inc1(:,ite,:) = Inc1R;
    R_inc2(:,ite,:) = Inc2R;
end
Scal1 = [U_scal1,B_scal1];
Scal2 = [U_scal2,B_scal2];

Inc1 = [L_inc1,R_inc1];
Inc2 = [L_inc2,R_inc2];


HH1 = Tensor_Mat_Mult(Pol_1,Inc1);
HH2 = Tensor_Mat_Mult(permute(Scal1,[2,1,3]),HH1);
Help1 = (mu_rel-1)*HH2;

II1 = Tensor_Mat_Mult(Pol_2,Inc2);
II2 = Tensor_Mat_Mult(permute(Scal2,[2,1,3]),II1);
Help2 = -kappa^2*(1-eps_rel)*II2;
integrand = (aa * bb)*pi*(Help1 + Help2);

% value of the integral

FarFieldMatrix = sum((wvec).*integrand,3);