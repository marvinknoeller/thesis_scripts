%> @file scatFromPCSphere.m
%> @brief Contains the scatFromPCSphere() function.
% ======================================================================
%> @brief For a given incident EntireWaveField as the electric field,
%> compute the scattered field for scattering from a perfectly conducting
%> sphere.
%>
%> @param[in] Einc Incident electric field, an EntireWaveField with the center 
%>    of the perfectly sphere as center of expansion.
%> @param[in] R Radius of the perfectly conducting sphere.
%> @param[out] Escat The scattered electric field, a RadiatingWaveField
function Escat = scatFromPCSphere(Einc, R)

alpha = zeros(1, (Einc.degree+1)^2 - 1 );
beta  = zeros(1, (Einc.degree+1)^2 - 1 );

for n=1:Einc.degree
    
    index = n^2:( (n+1)^2 - 1 );    
    
    j_n_over_kR = sqrt(pi / (2*Einc.wavenum) / R ) * besselj(n+1/2,Einc.wavenum*R) / (Einc.wavenum*R);
    j_n_min_1_R = sqrt(pi / (2*Einc.wavenum) ./ R ) .* besselj(n-1/2,Einc.wavenum*R);
    scale_M_nm = Einc.wavenum * R * j_n_over_kR;
    scale_curl_M_nm =  j_n_min_1_R - n * j_n_over_kR;
    
    h_n_over_kR = sqrt(pi / (2*Einc.wavenum) / R ) * besselh(n+1/2,1,Einc.wavenum*R) / (Einc.wavenum*R);
    h_n_min_1_R = sqrt(pi / (2*Einc.wavenum) ./ R ) .* besselh(n-1/2,1,Einc.wavenum*R);
    scale_N_nm = Einc.wavenum * R * h_n_over_kR;
    scale_curl_N_nm =  h_n_min_1_R - n * h_n_over_kR;
    
    alpha(index) = -Einc.alpha(index) * scale_M_nm / scale_N_nm;
    beta(index) = - Einc.beta(index) * scale_curl_M_nm / scale_curl_N_nm;
    
end

Escat = fields.RadiatingWaveField( Einc.wavenum, alpha, beta, Einc.z );

end