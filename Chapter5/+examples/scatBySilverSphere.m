function [ Einc, Escat, Etrans ] = scatBySilverSphere()

% set the parameters
lambda = 0.4;                % corresponds to 500 nm
R = 1;                       % sphere has a radius of 10 wavelengths
eps_r = -9.7999 + 0.31309i;  % value for silver at 500nm
                             % obtained from https://refractiveindex.info/
                             % with values from Johnson and Christy 1972
 
%eps_r = -2 + 0.05i;                             
                             
wavenum = 2*pi / lambda;

% set the incident field
d = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
p = [2; -1; -1];
z = [ 0; 0; 0 ];
PW = fields.PlaneWave( wavenum, d, p );

[ d_phi, d_th, ~ ] = cart2sph( d(1), d(2), d(3) );
d_th = pi/2 - d_th;

N = 25;
for n=1:N
    
    m_range = n^2:((n+1)^2 - 1);
    
    [U_d, V_d] = specialfunc.vectorSphericalHarmonics( d_th, d_phi, n );
    p_times_U_d = p(1) * U_d(:,:,1)' + p(2) * U_d(:,:,2)' + p(3) * U_d(:,:,3)';
    p_times_V_d = p(1) * V_d(:,:,1)' + p(2) * V_d(:,:,2)' + p(3) * V_d(:,:,3)';
    
    a_nm(m_range) = -4*pi * (1i)^n * p_times_V_d * exp( 1i * wavenum * dot(z,d) );
    b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d * exp( 1i * wavenum * dot(z,d) );
    
end

Einc = fields.EntireWaveField(wavenum, a_nm, b_nm, z);

% compute scattered and transmitted field
[Escat, Etrans] = scattering.scatFromPenetrableSphere(Einc, R, eps_r);

r = linspace(0,1,101);
phi = linspace(-pi,pi,201);
[RR, PHI] = meshgrid(r,phi);
X1 = RR .* cos(PHI);
X2 = RR .* sin(PHI);
X3 = zeros(size(X1));

figure(1);
[E1,E2,E3] = Etrans.eval(X1,X2,X3);
surf( X1,X2, real(E1) );
shading interp
view(2)
colorbar

r = linspace(1,2,151);
phi = linspace(-pi,pi,401);
[RR, PHI] = meshgrid(r,phi);
X1 = RR .* cos(PHI);
X2 = RR .* sin(PHI);
X3 = zeros(size(X1));

figure(2);
[E1,E2,E3] = Einc.eval(X1,X2,X3);
surf( X1,X2, real(E1) );
shading interp
view(2)
colorbar

[a, b] = Escat.farFieldCoeffs;
EFF = fields.L2tangField(a,b);
figure(3);
semilogy(1:(N+1)^2-1,abs(EFF.a),1:(N+1)^2-1,abs(EFF.b));

end