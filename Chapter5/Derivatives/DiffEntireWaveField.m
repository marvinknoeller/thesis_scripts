function [DMnm, DcurlMnm] = DiffEntireWaveField(kappa, n, m,  X1, X2, X3)

x1 = X1(:).';
x2 = X2(:).';
x3 = X3(:).';

[ phi, theta, r ] = cart2sph( x1, x2, x3 );
theta = pi/2 - theta;
phi = phi + 2*pi;
DMnm = zeros(3,length(X1),3);
DcurlMnm = zeros(3,length(X1),3);

rhat = [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];
thetahat = [cos(theta).*cos(phi); cos(theta).*sin(phi); -sin(theta)];
phihat = [-sin(phi); cos(phi); zeros(1,length(X1))];

% Initialize variables

absm = abs(m);
Cnm_tilde = sqrt( (2*n + 1)/(4*pi) * factorial(n-absm)/factorial(n+absm));
Cnm = 1/(sqrt(n*(n+1))) * Cnm_tilde;
Pn = zeros(absm+3,length(X1));
Pn(1:n+1,:) = legendre(n,cos(theta),'unnorm').*(-1).^(0:n)';

% for small arguments, we have to use asymptotic forms for the 
% Bessel functions
r_small = ( r < 1e-5 );
j_n_over_r(r_small) = sqrt(pi / 16 ) * (r(r_small)/2).^(n-1) * kappa^n / gamma(n+3/2);
j_n_over_r(~r_small)  = sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+1/2,kappa*r(~r_small)) ./ r(~r_small) ;

j_np1(~r_small)  = sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+3/2,kappa*r(~r_small)) ;
j_np1(r_small) = sqrt(pi / 16 ) /  gamma(n+5/2) * kappa^(n+1) * r(r_small).^(n+1) / (2^(n));

% for theta close to 0 or pi we have to use asymptotic forms of the
% functions which appear in the derivatives. For non-special theta, we
% consider theta_normal
theta_small = ( theta < 1e-10 );
theta_pi = ( abs(pi - theta) < 1e-10 );
theta_normal = ~(theta_small + theta_pi);
%----------------------------------------------------------------------------------------------------------------------------------------
%% Beginn Mnm
absm_cot_Pnm = zeros(1,length(X1));
absm_cot_Pnm(theta_small) = n*(n+1)/2 .* ( absm == 1 );
absm_cot_Pnm(theta_pi) = (-1)^(n) * n * (n + 1) / 2 .* ( absm ==1 ) ;
absm_cot_Pnm(theta_normal) = absm * cot(theta(theta_normal)) .* Pn(absm + 1,theta_normal);

i_m_cot_Pnm(theta_small) = 1i * m * n*(n+1)/2 .* ( absm == 1 );
i_m_cot_Pnm(theta_pi) = 1i * m * (-1)^n * n * (n + 1) / 2 .* ( absm ==1 ) ;
i_m_cot_Pnm(theta_normal) = 1i * m * cot(theta(theta_normal)) .* Pn(absm + 1, theta_normal);

absmterm(theta_small) = absm * 1/8 *(n+2) * (n+1) * n * (n-1) .* ( absm == 2 ) - n * (n+1)/2 .*(absm == 0) ;
absmterm(theta_pi) = (-1)^n * absm * 1/8 *(n+2) * (n+1) * n * (n-1) .* ( absm == 2 ) - (-1)^(n) * n * (n+1)/2 .*(absm == 0) ;
absmterm(theta_normal) = absm * (-1./(sin(theta(theta_normal)).^2) .* Pn(absm+1, theta_normal) + cot(theta(theta_normal))...
    .* (absm * cot(theta(theta_normal)) .* Pn(absm+1, theta_normal)-Pn(absm+2, theta_normal))) - ((absm+1) .* cot(theta(theta_normal))...
    .* Pn(absm+2,theta_normal) - Pn(absm+3,theta_normal)) ;

imterm(theta_small) = 1i * m * 1/8 *(n+2) * (n+1) * n * (n-1) .* ( absm == 2 ) ;
imterm(theta_pi) = 1i * m * (-1)^n * 1/8 *(n+2) * (n+1) * n * (n-1) .* (absm == 2 );
imterm(theta_normal) = 1i * m * (-1./(sin(theta(theta_normal)).^2) .* Pn(absm+1, theta_normal) + cot(theta(theta_normal))...
    .* (absm * cot(theta(theta_normal)) .* Pn(absm+1, theta_normal)-Pn(absm+2, theta_normal)));

phi_term1(theta_small) = -n * (n+1)/2 .* (1i * m * exp(1i*m*phi(theta_small)) .* sin(phi(theta_small)) + exp(1i*m*phi(theta_small)).*cos(phi(theta_small))) .* ( absm == 0 ) ...
    + (1i * m * 1/8 * (n+2) * (n+1) * n * (n-1) * exp(1i*m*phi(theta_small))...
    .* sin(phi(theta_small)) -2 * 1/8 * (n+2) * (n+1) * n * (n-1) *exp(1i*m*phi(theta_small)) .* cos(phi(theta_small))) .* ( absm == 2);
phi_term1(theta_pi) = -n * (n+1)/2 * (-1)^(n+1) .* (1i * m * exp(1i*m*phi(theta_pi)) .* sin(phi(theta_pi)) + exp(1i*m*phi(theta_pi)).*cos(phi(theta_pi))) .* ( absm == 0 ) ...
    + (1i * m * (-1)^(n+1) * 1/8 * (n+2) * (n+1) * n * (n-1) * exp(1i*m*phi(theta_pi)) .* sin(phi(theta_pi))...
    - 2 * (-1)^(n+1) * 1/8 * (n+2) * (n+1) * n * (n-1) *exp(1i*m*phi(theta_pi)) .* cos(phi(theta_pi))) .* ( absm == 2);
phi_term1(theta_normal) = 1./(sin(theta(theta_normal))) .* ( ( absm .* cot(theta(theta_normal)) .* Pn(absm+1, theta_normal) - Pn(absm+2, theta_normal))...
    .* (1i.*m.*exp(1i*m*phi(theta_normal)) .* sin(phi(theta_normal)) + exp(1i*m*phi(theta_normal)) .* cos(phi(theta_normal)))...
    + cot(theta(theta_normal)) .* 1i.* m .* Pn(absm +1,theta_normal)...
    .* (1i .* m .* exp(1i * m * phi(theta_normal)) .* cos(phi(theta_normal)) - exp(1i * m * phi(theta_normal)) .* sin(phi(theta_normal))));

phi_term2(theta_small) = n * (n+1)/2 .* (1i * m * exp(1i*m*phi(theta_small)) .* cos(phi(theta_small)) - exp(1i*m*phi(theta_small)).*sin(phi(theta_small))) .* ( absm == 0 )...
    - (1i * m * 1/8 * (n+2) * (n+1) * n * (n-1) * exp(1i*m*phi(theta_small))...
    .* cos(phi(theta_small)) -2 * 1/8 * (n+2) * (n+1) * n * (n-1) *exp(1i*m*phi(theta_small)) .* sin(phi(theta_small))) .* ( absm == 2);
phi_term2(theta_pi) = n * (n+1)/2 * (-1)^(n+1) .* (1i * m * exp(1i*m*phi(theta_pi)) .* cos(phi(theta_pi)) - exp(1i*m*phi(theta_pi)).*sin(phi(theta_pi))) .* ( absm == 0 )  ...
    - (1i * m * (-1)^(n+1) * 1/8 * (n+2) * (n+1) * n * (n-1) * exp(1i*m*phi(theta_pi)) .* cos(phi(theta_pi))...
    - 2 * (-1)^(n+1) * 1/8 * (n+2) * (n+1) * n * (n-1) *exp(1i*m*phi(theta_pi)) .* sin(phi(theta_pi))) .* ( absm == 2);
phi_term2(theta_normal) = 1./(sin(theta(theta_normal))) .* ( -( absm .* cot(theta(theta_normal)) .* Pn(absm+1, theta_normal) - Pn(absm+2, theta_normal))...
    .* (1i.*m.*exp(1i*m*phi(theta_normal)) .* cos(phi(theta_normal)) - exp(1i*m*phi(theta_normal)) .* sin(phi(theta_normal)))...
    +cot(theta(theta_normal)) .* 1i.* m .* Pn(absm +1,theta_normal)...
    .* (1i .* m .* exp(1i * m * phi(theta_normal)) .* sin(phi(theta_normal)) + exp(1i * m * phi(theta_normal)) .* cos(phi(theta_normal))));

phi_term3(theta_small) = m^2 * n * (n+1) /2 .* ( absm == 1);
phi_term3(theta_pi) = m^2 * (-1)^(n+1) * n * (n+1) /2 .* (absm ==1);
phi_term3(theta_normal) = m^2 .* 1./sin(theta(theta_normal)) .* Pn( absm +1 , theta_normal);

%% Ableitung nach r
dmn1_r = kappa * (-j_np1 + n/(kappa ) .* j_n_over_r ) .* Cnm .*( ( absm_cot_Pnm - Pn(absm + 2,:) ) .* exp(1i*m*phi) .* sin(phi) ...
    + i_m_cot_Pnm .* exp(1i*m*phi) .* cos(phi)) .* rhat;
dmn2_r = kappa * (-j_np1 + n/(kappa ) .* j_n_over_r ) .* Cnm .*(-( absm_cot_Pnm - Pn(absm + 2,:) ) .* exp(1i*m*phi) .* cos(phi) ...
    + i_m_cot_Pnm .* exp(1i*m*phi) .* sin(phi)) .* rhat;
dmn3_r = -kappa * (-j_np1 + n/(kappa ) .* j_n_over_r ) .* Cnm .* 1i .* m .* Pn(absm + 1,:) .* exp(1i*m*phi) .* rhat;

%% Ableitung nach phi
dmn1_phi = j_n_over_r .* Cnm .* phi_term1 .* phihat;
dmn2_phi = j_n_over_r .* Cnm .* phi_term2 .* phihat;
dmn3_phi = j_n_over_r .* Cnm .* phi_term3 .* exp(1i*m*phi) .* phihat;

%% Ableitung nach theta
dmn1_theta = j_n_over_r .* Cnm .* ( absmterm .* exp(1i*m*phi) .* sin(phi) ...
    + imterm .* exp(1i*m*phi) .* cos(phi)) .* thetahat;
dmn2_theta = j_n_over_r .* Cnm .* (-absmterm .* exp(1i*m*phi) .* cos(phi) ...
    + imterm .* exp(1i*m*phi) .* sin(phi)) .* thetahat;
dmn3_theta = -j_n_over_r .* Cnm .* 1i .* m .* (absm_cot_Pnm - Pn(absm+2,:)) .* exp(1i*m*phi) .* thetahat; %% Zusammen setzen und testen

%% Zusammensetzen
DMnm1 = dmn1_r + dmn1_phi + dmn1_theta;
DMnm2 = dmn2_r + dmn2_phi + dmn2_theta;
DMnm3 = dmn3_r + dmn3_phi + dmn3_theta;

DMnm(1,:,:) = DMnm1.';
DMnm(2,:,:) = DMnm2.';
DMnm(3,:,:) = DMnm3.';

%----------------------------------------------------------------------------------------------------------------------------------------
%% Beginn CurlMnm
Pndd = (n+2) * (n+1) * n * (n-1) / 8;

% only for n>1 (n-1)j_n/r^2
nm1_j_n_over_r2(r_small) = (n-1) * sqrt(pi / 16 ) * (1/2).^(n-1) * r(r_small).^(n-2) * kappa^n / gamma(n+3/2);
nm1_j_n_over_r2(~r_small)  = (n-1) * sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+1/2,kappa*r(~r_small)) ./ (r(~r_small)).^2 ;
if n==1
    nm1_j_n_over_r2(r_small) = 0;
    nm1_j_n_over_r2(~r_small)  = 0;
end
   
% It's a trick : jn_over_r2 is singular, but the terms in the derivative
% w.r.t. phi cancel out - so I can set it to 0 
Pj_n_over_r2(r_small) = sqrt(pi / 16 ) * (1/2).^(n-1) * r(r_small).^(n-2) * kappa^n / gamma(n+3/2);
Pj_n_over_r2(~r_small)  = sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+1/2,kappa*r(~r_small)) ./ (r(~r_small)).^2 ;
if n==1
    Pj_n_over_r2(r_small) = 0;
    Pj_n_over_r2(~r_small)  = 0;
end

j_np2(~r_small)  = sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+5/2,kappa*r(~r_small)) ;
j_np2(r_small) = sqrt(pi / 16 ) /  gamma(n+7/2) * kappa^(n+2) * r(r_small).^(n+2) / (2^(n+1));
j_np1_over_r(~r_small)  = sqrt(pi / (2*kappa) ./ r(~r_small) ) .* besselj(n+3/2,kappa*r(~r_small))./ r(~r_small) ;
j_np1_over_r(r_small) = sqrt(pi / 16 ) /  gamma(n+5/2) * kappa^(n+1) * r(r_small).^(n) / (2^(n));

r_term1 = -kappa * j_np1_over_r + nm1_j_n_over_r2;
r_term2 = (n+1) * r_term1 - kappa^2 * (-j_np2 + (n+1)/kappa * j_np1_over_r);
r_term3 = r_term2; 

theta_term1 = Pn(absm + 1, :) .* sin(theta);
theta_term2(theta_small) = n*(n+1)/2 .* ( absm == 1 );
theta_term2(theta_pi) = (-1).^(n+1) .* n*(n+1)/2 .* ( absm == 1 );
theta_term2(theta_normal) = absm * cot(theta(theta_normal)) .* Pn(absm + 1, theta_normal) .* cos(theta(theta_normal)) - Pn(absm +2, theta_normal) .* cos(theta(theta_normal));
theta_term3(theta_small) = 1i * m * n*(n+1)/2 .* ( absm == 1 );
theta_term3(theta_pi) = 1i * m * (-1)^(n+1) * n*(n+1)/2 .* ( absm == 1 );
theta_term3(theta_normal) = 1i * m * 1./sin(theta(theta_normal)) .* Pn(absm + 1, theta_normal); 
theta_term4 = Pn(absm + 1, :) .* cos(theta);
theta_term5(theta_normal) = absm * cot(theta(theta_normal)) .* Pn(absm + 1, theta_normal) .* (-sin(theta(theta_normal))) - Pn(absm +2, theta_normal) .* (-sin(theta(theta_normal)));
theta_term5(~theta_normal) = 0;

absmm1(theta_small) = (absm -1) * Pndd .* (absm == 2);
absmm1(theta_pi) = (absm -1) * (-1)^n * Pndd .* (absm == 2);
absmm1(theta_normal) = (absm-1) * 1./sin(theta(theta_normal)).^2 .* Pn(absm +1, theta_normal);

absmmm2(theta_small) = (absm -m^2) * Pndd .* (absm == 2);
absmmm2(theta_pi) = (absm -m^2) * (-1)^n * Pndd .* (absm == 2);
absmmm2(theta_normal) = (absm-m^2) * 1./sin(theta(theta_normal)).^2 .* Pn(absm +1, theta_normal);

cotPnmp1(theta_small) = n * (n+1) / 2 .* ( absm == 0 );
cotPnmp1(theta_pi) = (-1)^n * n * (n+1) / 2 .* ( absm == 0 );
cotPnmp1(theta_normal) = cot(theta(theta_normal)) .* Pn(absm + 2, theta_normal) ; %%%%%%%%%%

phi_term1 = Pj_n_over_r2 .* ( (n - absm) * Pn(absm + 1,:) + absmm1 - cotPnmp1) ; 
phi_term2 = - kappa * j_np1_over_r .* (absmm1 - absm * Pn(absm + 1, :) - cotPnmp1);
phi_term3 = Pj_n_over_r2 .* ( (n - absm) * Pn(absm + 1,:) + absmmm2 - cotPnmp1);
phi_term4 = - kappa * j_np1_over_r .* (absmmm2 - absm * Pn(absm + 1, :) - cotPnmp1);

absm_cot_Pnm(theta_small) = n*(n+1)/2 .* ( absm == 1 );
absm_cot_Pnm(theta_pi) = (-1)^(n) * n * (n + 1) / 2 .* ( absm ==1 ) ;
absm_cot_Pnm(theta_normal) = absm * cot(theta(theta_normal)) .* Pn(absm + 1,theta_normal);

nmabsm_cot_Pnm(theta_small) = ( n - absm) * n*(n+1)/2 .* ( absm == 1 );
nmabsm_cot_Pnm(theta_pi) = (n - absm) * (-1)^n * n*(n+1)/2 .* ( absm == 1);
nmabsm_cot_Pnm(theta_normal) = (n - absm) * cot(theta(theta_normal)) .* Pn(absm + 1, theta_normal);

phi_term5 = Pj_n_over_r2 .* ( nmabsm_cot_Pnm + Pn(absm + 2,:));
phi_term6 = - kappa * j_np1_over_r .* (-absm_cot_Pnm + Pn(absm + 2,:));

absm_over_sin2_term(theta_small) = absm * Pndd .* ( absm ==2 );
absm_over_sin2_term(theta_pi) = absm * (-1)^(n+1) * Pndd .* ( absm ==2 );
absm_over_sin2_term(theta_normal) = absm./sin(theta(theta_normal)).^2 .* cos(theta(theta_normal)) .* Pn( absm + 1, theta_normal);

absm2_cot2_term(theta_small) = absm^2 * Pndd .* ( absm ==2 );
absm2_cot2_term(theta_pi) = absm^2 * (-1)^(n+1) * Pndd .* ( absm ==2 );
absm2_cot2_term(theta_normal) = absm^2 .* cot(theta(theta_normal)).^2 .* Pn(absm+1 , theta_normal) .* cos(theta(theta_normal));

absmp1_cot_term(theta_small) = ( absm + 1 ) .* n*(n+1)/2 .* ( absm ==0 );
absmp1_cot_term(theta_pi) = ( absm + 1 ) .* (-1)^(n+1) .* n*(n+1)/2 .* ( absm ==0 );
absmp1_cot_term(theta_normal) = ( absm + 1 ) .* cot(theta(theta_normal)) .* cos(theta(theta_normal)) .* Pn( absm + 2, theta_normal);

over_sin2_term(theta_small) =  Pndd .* ( absm ==2 );
over_sin2_term(theta_pi) =  (-1)^(n+1) * Pndd .* ( absm ==2 );
over_sin2_term(theta_normal) = 1./sin(theta(theta_normal)).^2 .* cos(theta(theta_normal)) .* Pn( absm + 1, theta_normal);

over_sin_term(theta_small) = n*(n+1)/2 .* ( absm == 0 );
over_sin_term(theta_pi) = (-1)^(n+1) .* n*(n+1)/2 .* ( absm == 0 );
over_sin_term(theta_normal) = 1./sin(theta(theta_normal)) .* Pn( absm + 2, theta_normal);

absm_cot_Pnmp1_cos(theta_small) = 0;
absm_cot_Pnmp1_cos(theta_pi) = 0;
absm_cot_Pnmp1_cos(theta_normal) = absm .* Pn( absm + 2,theta_normal) .* cos(theta(theta_normal)) .* cot(theta(theta_normal));

th_term1 = n * (absm + 1) .* cos(theta) .* Pn( absm + 1, :) - n * sin(theta) .* Pn( absm + 2, :);
th_term2 = -absm_over_sin2_term + absm2_cot2_term...
    - absm_cot_Pnmp1_cos... 
    - absmp1_cot_term + Pn( absm + 3, :) .* cos(theta) ...
    - absm_cot_Pnm.*sin(theta) + Pn( absm + 2, :) .* sin(theta);
th_term3 = ( absm - 1 ) .* over_sin2_term - over_sin_term;

over_sin_term2(theta_small) = n*(n+1)/2 .* ( absm == 1 );
over_sin_term2(theta_pi) = (-1)^(n+1) .* n*(n+1)/2 .* ( absm == 1 );
over_sin_term2(theta_normal) = 1./sin(theta(theta_normal)) .* Pn( absm + 1, theta_normal); 

th_term4 = n * ( absm_cot_Pnm - Pn( absm + 2, :)) .* cos(theta) - n * Pn( absm + 1, :) .* sin(theta);
th_term5 = absm * ( over_sin_term2 - cos(theta) .* ( absm_cot_Pnm - Pn( absm +2, :) ))...
    + ( (absm + 1) * cos(theta) .* Pn( absm + 2, :) - sin(theta) .* Pn( absm +3, :))...
    -( cos(theta) .* absm_cot_Pnm - Pn( absm + 2, :) .* cos(theta));
%% Ableitung nach r
dmn1_curl_r = (sqrt(n * (n+1)) * Cnm_tilde .* (r_term1 .* theta_term1 .* exp(1i * m * phi) .* cos(phi)) ...
    + Cnm .* r_term2 .* theta_term2 .* exp(1i * m * phi) .* cos(phi) ...
    + Cnm .* r_term3 .* theta_term3 .* exp(1i * m * phi) .* (-sin(phi)) )  .* rhat;
dmn2_curl_r = (sqrt(n * (n+1)) * Cnm_tilde .* (r_term1 .* theta_term1 .* exp(1i * m * phi) .* sin(phi)) ...
    + Cnm .* r_term2 .* theta_term2 .* exp(1i * m * phi) .* sin(phi) ...
    + Cnm .* r_term3 .* theta_term3 .* exp(1i * m * phi) .* (cos(phi)) )  .* rhat;
dmn3_curl_r = (sqrt(n * (n+1)) * Cnm_tilde .* (r_term1 .* theta_term4 .* exp(1i * m * phi) ) ...
    + Cnm .* r_term2 .* theta_term5 .* exp(1i * m * phi) ) .* rhat;

%% Ableitung nach phi
dmn1_curl_phi = ((Cnm * (n+1) * 1i * m * phi_term1 + Cnm * 1i * m * phi_term2) .* exp(1i * m * phi) .* cos(phi)...
    + ( -Cnm * (n+1) * phi_term3 - Cnm * phi_term4 ) .* exp(1i * m * phi) .* sin(phi) ) .* phihat;
dmn2_curl_phi = ((Cnm * (n+1) * 1i * m * phi_term1 + Cnm * 1i * m * phi_term2) .* exp(1i * m * phi) .* sin(phi)...
    + ( -Cnm * (n+1) * phi_term3 - Cnm * phi_term4 ) .* exp(1i * m * phi) .* (-cos(phi)) ) .* phihat;
dmn3_curl_phi = (Cnm * 1i * m * (n+1) * phi_term5 .* exp(1i * m * phi)... 
    + Cnm * 1i * m * phi_term6 .* exp(1i * m * phi)) .* phihat;
%% Ableitung nach theta
dmn1_curl_theta = (Cnm * (n+1) * Pj_n_over_r2 .* (th_term1 + th_term2) .* exp(1i * m * phi) .* cos(phi)...
    + Cnm * (- kappa * j_np1_over_r) .* th_term2 .* exp(1i * m * phi) .* cos(phi)...
    + Cnm * ( (n+1)* Pj_n_over_r2 - kappa * j_np1_over_r ) .* (1i * m * th_term3) .* exp(1i * m * phi) .* (-sin(phi)) ) .*thetahat;
dmn2_curl_theta = (Cnm * (n+1) * Pj_n_over_r2 .* (th_term1 + th_term2) .* exp(1i * m * phi) .* sin(phi)...
    + Cnm * (- kappa * j_np1_over_r) .* th_term2 .* exp(1i * m * phi) .* sin(phi)...
    + Cnm * ( (n+1)* Pj_n_over_r2 - kappa * j_np1_over_r ) .* (1i * m * th_term3) .* exp(1i * m * phi) .* (cos(phi)) ) .*thetahat;
dmn3_curl_theta = (Cnm * (n+1) * Pj_n_over_r2 .* (th_term4 + th_term5) .* exp(1i * m * phi)...
    + Cnm * (- kappa * j_np1_over_r) .* th_term5 .* exp(1i * m * phi)) .* thetahat;
%% Zusammensetzen

DcurlMnm1 = dmn1_curl_r + dmn1_curl_phi + dmn1_curl_theta;
DcurlMnm2 = dmn2_curl_r + dmn2_curl_phi + dmn2_curl_theta;
DcurlMnm3 = dmn3_curl_r + dmn3_curl_phi + dmn3_curl_theta;

DcurlMnm(1,:,:) = DcurlMnm1.';
DcurlMnm(2,:,:) = DcurlMnm2.';
DcurlMnm(3,:,:) = DcurlMnm3.';


end

