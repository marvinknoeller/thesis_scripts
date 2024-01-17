function der_ffm = derive_farfieldmatrix_SplineRotationRMFFD(p,h,Var,R,S,T)

eps_rel = Var.eps_rel;
mu_rel = Var.mu_rel;
kappa = Var.kappa;
M = Var.M;
num_x = size(p,2);
N = Var.N;
for nn=1:N
    nvec(nn^2:nn^2+2*nn)=nn;
end
% size of vectors and matrix
Q = 2*N*(N+2);
Qd2 = Q/2;

aa = Var.aa;
bb = Var.bb;



const_eps1 = (aa + bb) / (aa + eps_rel*bb);
const_eps2 = (aa + bb) / (bb + eps_rel*aa);

const_mu1 = (aa + bb) / (aa + mu_rel*bb);
const_mu2 = (aa + bb) / (bb + mu_rel*aa);

n = size(p,2) - 1;
num_x = n+1;
lenV = n*M-(n-1);

%%
L_inc1 = zeros(3,Qd2,lenV);
L_inc2 = zeros(3,Qd2,lenV);
R_inc1 = zeros(3,Qd2,lenV);
R_inc2 = zeros(3,Qd2,lenV);
U_scal1 = zeros(3,Qd2,lenV);
U_scal2 = zeros(3,Qd2,lenV);
B_scal1 = zeros(3,Qd2,lenV);
B_scal2 = zeros(3,Qd2,lenV);
DiffU_scal1 = zeros(3,Qd2,lenV);
DiffU_scal2 = zeros(3,Qd2,lenV);
DiffB_scal1 = zeros(3,Qd2,lenV);
DiffB_scal2 = zeros(3,Qd2,lenV);
DiffL_inc1 = zeros(3,Qd2,lenV);
DiffL_inc2 = zeros(3,Qd2,lenV);
DiffR_inc1 = zeros(3,Qd2,lenV);
DiffR_inc2 = zeros(3,Qd2,lenV);

N = -1 + sqrt( 1 + Q/2 );
nvec = [];
mvec = [];
for ki = 1 : N
    nvec = [nvec, ki*ones(1,2*ki+1)];
    mvec = [mvec, [-ki:1:ki]];
end
D_mu = diag([1,const_mu1, const_mu2]);

[~,~,coefs,~,ts] = splinepoints(p,Var.M);
%%%%%%%%%%%
[p_in_between,der_p,derder_p,tt] = allpoints(coefs,ts,num_x,M);
delta_t = tt(2:end) - tt(1:end-1);

[~,~,coefsh,~,tsh] = splinepoints(h,Var.M,ts);
[h_in_between,der_h,derder_h,~] = allpoints(coefsh,tsh,num_x,M);

cc = p + h;
[~,~,coefscc,~,tsh] = splinepoints(cc,Var.M,ts);
[cc_in_between,der_cc,derder_h,~] = allpoints(coefscc,tsh,num_x,M);
for ell = 1 : length(p_in_between)
    Tnew(:,ell) = der_cc(:,ell)/norm(der_cc(:,ell));
    Rnew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * R(:,ell) ...
        - dotReal(S(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
        - dotReal(R(:,ell),Tnew(:,ell)) * T(:,ell);
    Snew(:,ell) = dotReal(T(:,ell),Tnew(:,ell)) * S(:,ell) ...
        + dotReal(R(:,ell),Tnew(:,ell)) / (1+ dotReal(T(:,ell),Tnew(:,ell))) * cross(T(:,ell),Tnew(:,ell)) ...
        - dotReal(S(:,ell),Tnew(:,ell)) * T(:,ell);
end

for ii = 1:length(tt)
        DDreal(:,:,ii) = [real(const_eps1), 0; 0 real(const_eps2)];
        DDimag(:,:,ii) = [imag(const_eps1), 0; 0 imag(const_eps2)];
end


%%
ww = zeros(1,size(p_in_between,2));
tanvec = zeros(3,size(p_in_between,2));
lenvec = zeros(1,size(p_in_between,2));
Pol_1 = zeros(3,3,size(p_in_between,2));
Pol_2 = zeros(3,3,size(p_in_between,2));
der_Pol_1 = zeros(3,3,size(p_in_between,2));
der_Pol_2 = zeros(3,3,size(p_in_between,2));
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
    lenvec(ki) = norm(der_p(:,ki));
    %the matrix V, we need to diag. the pol. tensor
    tan = T(:,ki);
    tanvec(:,ki) = tan;
    nor = R(:,ki);
    bnor = S(:,ki);
    Vmu = [tan, nor, bnor];
    VVV = [tan, nor, bnor];
    WWW = [tan, nor, bnor];

    Pol_1(:,:,ki) = Vmu*diag([1,const_mu1,const_mu2])*(Vmu.');
    Pol_2(:,:,ki) = VVV*diag([1,DDreal(1,1,ki),DDreal(2,2,ki)])*(VVV.') ...
        +1i*WWW*diag([0,DDimag(1,1,ki),DDimag(2,2,ki)])*(WWW.');
    Tdiffn(:,ki) = 1/(norm(der_p(:,ki)))*(dotReal(der_h(:,ki),nor)*nor+dotReal(der_h(:,ki),bnor)*bnor);
    Rdiffn(:,ki) = - 1/(norm(der_p(:,ki))) * dotReal(der_h(:,ki),nor)*tan;
    Sdiffn(:,ki) = - 1/(norm(der_p(:,ki))) * dotReal(der_h(:,ki),bnor)*tan;
    nTbT = [Rdiffn(:,ki), Sdiffn(:,ki)];
    V_ph = [Tdiffn(:,ki)...
        , nTbT];
    der_Pol_1(:,:,ki) = (V_ph*D_mu*(Vmu.')+Vmu*D_mu*(V_ph.'));
    der_Pol_2(:,:,ki) = (V_ph*diag([1,DDreal(1,1,ki),DDreal(2,2,ki)])*(VVV.')...
        +VVV*diag([1,DDreal(1,1,ki),DDreal(2,2,ki)])*(V_ph.'))...
        +1i * (V_ph*diag([0,DDimag(1,1,ki),DDimag(2,2,ki)])*(WWW.')...
        +WWW*diag([0,DDimag(1,1,ki),DDimag(2,2,ki)])*(V_ph.'));
end
wvec = permute(repmat(ww,Q,1,Q),[1,3,2]); %Q
X1 = p_in_between(1,:);
X2 = p_in_between(2,:);
X3 = p_in_between(3,:);


% Consider terms der_ffm_AB at position (ite1 , ite2)
for ite = 1: Qd2
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
    [sMnmU,scurlMnmU] = DiffEntireWaveField(kappa,nvec(ite),mvec(ite),X1,X2,X3); %scalar product Mnm and scalar product curlMnm

    persMnmU = permute(sMnmU,[1,3,2]);
    perscurlMnmU = permute(scurlMnmU,[1,3,2]);
    U_DMnm_scal1 = zeros(3,lenV);
    U_curlDMnm_scal2 = zeros(3,lenV);
    for kk = 1 : size(h_in_between,2)
        U_DMnm_scal1(:,kk) = persMnmU(:,:,kk) * h_in_between(:,kk);%%%%
        U_curlDMnm_scal2(:,kk) = perscurlMnmU(:,:,kk) * h_in_between(:,kk);
    end
    U_DMnm_scal1 = 4*pi*kappa*1/((1i)^(n_tilde-1))*conj(U_DMnm_scal1);
    U_curlDMnm_scal2 = (4*pi/kappa)*1/((1i)^(n_tilde-1))*conj(U_curlDMnm_scal2);
    DiffU_scal1(:,ite,:) = U_DMnm_scal1;
    DiffU_scal2(:,ite,:) = U_curlDMnm_scal2;
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
    [sMnmB,scurlMnmB] = DiffEntireWaveField(kappa,nvec(ite),mvec(ite),X1,X2,X3); %scalar product Mnm and scalar product curlMnm

    persMnmB = permute(sMnmB,[1,3,2]);
    perscurlMnmB = permute(scurlMnmB,[1,3,2]);
    B_DMnm_scal1 = zeros(3,lenV);
    B_curlDMnm_scal2 = zeros(3,lenV);
    for kk = 1 : size(h_in_between,2)
        B_DMnm_scal1(:,kk) = persMnmB(:,:,kk) * h_in_between(:,kk);
        B_curlDMnm_scal2(:,kk) = perscurlMnmB(:,:,kk) * h_in_between(:,kk);
    end
    B_DMnm_scal1 = -4*pi/((1i)^(n_tilde))*conj(B_DMnm_scal1);
    B_curlDMnm_scal2 = -4*pi/((1i)^(n_tilde))*conj(B_curlDMnm_scal2);
    DiffB_scal1(:,ite,:) = B_curlDMnm_scal2;
    DiffB_scal2(:,ite,:) = B_DMnm_scal1;
    
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
    [iMnmL,icurlMnmL] = DiffEntireWaveField(kappa,nvec(ite),mvec(ite),X1,X2,X3); %scalar product Mnm and scalar product curlMnm

    persMnmL = permute(iMnmL,[1,3,2]);
    perscurlMnmL = permute(icurlMnmL,[1,3,2]); %[1,3,2]
    L_DMnm_inc1 = zeros(3,lenV);
    L_curlDMnm_inc2 = zeros(3,lenV);
    for kk = 1 : size(h_in_between,2)
        L_DMnm_inc1(:,kk) = persMnmL(:,:,kk) * h_in_between(:,kk);
        L_curlDMnm_inc2(:,kk) = perscurlMnmL(:,:,kk) * h_in_between(:,kk);
    end
    L_DMnm_inc1 = -4*pi*kappa*((1i)^(n_tilde+1))*(L_DMnm_inc1);
    L_curlDMnm_inc2 = (4*pi/kappa)*((1i)^(n_tilde-1))*(L_curlDMnm_inc2);
    DiffL_inc1(:,ite,:) = L_DMnm_inc1;
    DiffL_inc2(:,ite,:) = L_curlDMnm_inc2;
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
    [iMnmR,icurlMnmR] = DiffEntireWaveField(kappa,nvec(ite),mvec(ite),X1,X2,X3); %scalar product Mnm and scalar product curlMnm
    persMnmR = permute(iMnmR,[1,3,2]);
    perscurlMnmR = permute(icurlMnmR,[1,3,2]);
    R_DMnm_inc1 = zeros(3,lenV);
    R_curlDMnm_inc2 = zeros(3,lenV);
    for kk = 1 : size(h_in_between,2)
        R_DMnm_inc1(:,kk) = persMnmR(:,:,kk) * h_in_between(:,kk);
        R_curlDMnm_inc2(:,kk) = perscurlMnmR(:,:,kk) * h_in_between(:,kk);
    end
    R_DMnm_inc1 = -4*pi*((1i)^(n_tilde))*(R_DMnm_inc1);
    R_curlDMnm_inc2 = -4*pi*((1i)^(n_tilde))*(R_curlDMnm_inc2);
    DiffR_inc1(:,ite,:) = R_curlDMnm_inc2;
    DiffR_inc2(:,ite,:) = R_DMnm_inc1;
end
Scal1 = [U_scal1,B_scal1];
Scal2 = [U_scal2,B_scal2];
DiffScal1 = [DiffU_scal1,DiffB_scal1];
DiffScal2 = [DiffU_scal2,DiffB_scal2];

Inc1 = [L_inc1,R_inc1];
Inc2 = [L_inc2,R_inc2];
DiffInc1 = [DiffL_inc1,DiffR_inc1];
DiffInc2 = [DiffL_inc2,DiffR_inc2];

%% integrand1
HH1 = Tensor_Mat_Mult(Pol_1,Inc1);
HH2 = Tensor_Mat_Mult(permute(DiffScal1,[2,1,3]),HH1);
integrand1 = HH2;
%% integrand2
HH1 = Tensor_Mat_Mult(der_Pol_1,Inc1);
HH2 = Tensor_Mat_Mult(permute(Scal1,[2,1,3]),HH1);
integrand2 = HH2;
%% integrand3
HH1 = Tensor_Mat_Mult(Pol_1,DiffInc1);
HH2 = Tensor_Mat_Mult(permute(Scal1,[2,1,3]),HH1);
integrand3 = HH2;
%% integrand4
jj = (dotReal(tanvec,der_h)./lenvec);
jjmat = permute(repmat(jj,Q,1,Q),[1,3,2]);
HH1 = Tensor_Mat_Mult(Pol_1,Inc1);
HH2 = Tensor_Mat_Mult(permute(Scal1,[2,1,3]),HH1) .* jjmat;
integrand4 = HH2;
%%
Help1 = aa*bb * pi * (mu_rel-1) * (integrand1 + integrand2 + integrand3 + integrand4);

%% integrand5
HH1 = Tensor_Mat_Mult(Pol_2,Inc2);
HH2 = Tensor_Mat_Mult(permute(DiffScal2,[2,1,3]),HH1);
integrand5 = HH2;
%% integrand6
HH1 = Tensor_Mat_Mult(der_Pol_2,Inc2);
HH2 = Tensor_Mat_Mult(permute(Scal2,[2,1,3]),HH1);
integrand6 = HH2;
%% integrand7
HH1 = Tensor_Mat_Mult(Pol_2,DiffInc2);
HH2 = Tensor_Mat_Mult(permute(Scal2,[2,1,3]),HH1); 
integrand7 = HH2;
%% integrand8
jj = (dotReal(tanvec,der_h)./lenvec);
jjmat = permute(repmat(jj,Q,1,Q),[1,3,2]);
HH1 = Tensor_Mat_Mult(Pol_2,Inc2);
HH2 = Tensor_Mat_Mult(permute(Scal2,[2,1,3]),HH1) .* jjmat;
integrand8 = HH2;
%%
Help2 = -aa * bb * pi * kappa^2 * (1-eps_rel) * (integrand5 + integrand6 + integrand7 + integrand8);

der_ffm = sum((wvec) .* ( Help1 + Help2),3);
end