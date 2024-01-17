%> @file scatFromPenetrableSphere.m
%> @brief Contains the scatFromPenetrableSphere() function.
% ======================================================================
%> @brief For a given EntireWaveField as the incident electric field,
%> compute the scattered field for scattering from a penetrable homogenous 
%> sphere.
%>
%> @param[in] Einc Incident electric field, an EntireWaveField with the center 
%>    of the sphere as center of expansion.
%> @param[in] R Radius of the penetrable sphere.
%> @param[in] eps_r relative permittivity inside the sphere.
%> @param[out] Escat The scattered electric field, a RadiatingWaveField
%> @param[out] Etrans The transmitted electric field, an EntireWaveField
function [Escat, Etrans] = scatFromPenetrableSphere(Einc, R, eps_r)

alpha = zeros(1, (Einc.degree+1)^2 - 1 );
beta  = zeros(1, (Einc.degree+1)^2 - 1 );
gamma  = zeros(1, (Einc.degree+1)^2 - 1 );
delta  = zeros(1, (Einc.degree+1)^2 - 1 );

wavenum_inside = sqrt(eps_r) * Einc.wavenum;

for n=1:Einc.degree
    
    index = n^2:( (n+1)^2 - 1 );    
    
    j_n_out = sqrt(pi / (2*Einc.wavenum) / R ) * besselj(n+1/2,Einc.wavenum*R);
    j_n_min_1_out = sqrt(pi / (2*Einc.wavenum) ./ R ) .* besselj(n-1/2,Einc.wavenum*R);
    J_n_out = Einc.wavenum * R .* j_n_min_1_out - n * j_n_out;
    
    h_n = sqrt(pi / (2*Einc.wavenum) / R ) * besselh(n+1/2,1,Einc.wavenum*R);
    h_n_min_1 = sqrt(pi / (2*Einc.wavenum) ./ R ) .* besselh(n-1/2,1,Einc.wavenum*R);
    H_n = Einc.wavenum * R .* h_n_min_1 - n * h_n;
    
    j_n_in=  sqrt(pi / (2*wavenum_inside) / R ) * besselj(n+1/2,wavenum_inside*R);
    j_n_min_1_in = sqrt(pi / (2*wavenum_inside) ./ R ) .* besselj(n-1/2,wavenum_inside*R);
    J_n_in = wavenum_inside * R .* j_n_min_1_in - n * j_n_in;
    
    D1_n = H_n .* j_n_in - h_n .* J_n_in;
    D2_n = H_n .* j_n_in - 1/eps_r * h_n .* J_n_in;
    
    %alpha(index) = -Einc.alpha(index) * scale_M_nm / scale_N_nm;
    %beta(index) = - Einc.beta(index) * scale_curl_M_nm / scale_curl_N_nm;
            
    alpha(index) = Einc.alpha(index) .* (j_n_out .* J_n_in - J_n_out .* j_n_in) ./ D1_n;
    beta(index) = Einc.beta(index) .* ( 1/eps_r * j_n_out .* J_n_in - J_n_out .* j_n_in) ./ D2_n;
    
    gamma(index) = Einc.alpha(index) .* (j_n_out .* H_n - J_n_out .* h_n) ./ D1_n;
    delta(index) = Einc.beta(index) .* ( j_n_out .* H_n - J_n_out .* h_n) ./ D2_n / sqrt(eps_r);
    
end

Escat = fields.RadiatingWaveField( Einc.wavenum, alpha, beta, Einc.z );
Etrans = fields.EntireWaveField( wavenum_inside, gamma, delta, Einc.z );

end