function ge = FarField_Pert_Maxwell_E_spline(Var,p_in_between,Pol_1,Pol_2,ww,Z)
%% store the given parameters in local variables
%radius
roh=Var.roh;
%coefficients
eps_rel = Var.eps_rel;
mu_rel = Var.mu_rel;
eps0 = Var.eps0;
%direction of incidence
theta=Var.theta;
%compute the wave number
KA=Var.KA;
%polarization
A = Var.A;
%% Some help variables
% ge is the vector that contains the values of the pert. formula at the
% 2N(N-1) measure points. Note that it is complex valued.
ge = zeros(size(Z,2),3); 
w = ww;
%copy theta depending on the number of integration points
theta=theta*ones(1,length(p_in_between));
%incident wave
E0_f = @(X)sqrt(eps0)*A*exp(1i*KA*dotReal(X,theta)); 
E0(:,1,:) = E0_f(p_in_between);
%curl of incident wave
curl_E0_f = @(X)sqrt(eps0)*1i*KA*(cross(theta,A*ones(1,length(p_in_between)))).*exp(1i*KA*dotReal(X,theta));
curl_E0(:,1,:) = curl_E0_f(p_in_between); 
%% Compute the perturbation formula
% loop over the measurement points
for j=1:size(Z,2)
    x_hat = Z(:,j);
    copycross1 = repmat(cross([x_hat x_hat x_hat],eye(3)),1,1,length(p_in_between));
    Mat1 = Tensor_Mat_Mult(copycross1,Pol_1);
    copycross2 = repmat(cross([x_hat x_hat x_hat],(cross(eye(3),[x_hat x_hat x_hat]))),1,1,length(p_in_between));
    Mat2 = Tensor_Mat_Mult(copycross2,Pol_2);
    x_hat = x_hat*ones(1,length(p_in_between));
    % first term in the first integral
    fact = (exp(-1i*KA*dotReal(x_hat,p_in_between)));
    %1st matrix depends on measure point
    hh1=[fact;fact;fact].*squeeze((Tensor_Mat_Mult(Mat1,curl_E0)));
    %first integrand
    Help1 = 1i*KA*(mu_rel-1)*hh1;
    %second integrand
    Help2 = -KA^2*(1-eps_rel)*([fact;fact;fact].*squeeze((Tensor_Mat_Mult(Mat2,E0))));       
    % the integrand at the mesh points
    integrand = (roh^2)*pi*(Help1 + Help2);
    % value of the integral
    ge(j,:) = sum([w;w;w].*integrand,2);
end