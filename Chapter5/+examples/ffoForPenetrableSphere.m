%> @file ffoForPenetrableSphere.m
%> @brief Contains the ffoForPenetrableSphere() function.
% ======================================================================
%> @brief For a given penetrabe sphere, the far field operator is computed.
%>
%> @param[in] N maximal degree of vector spherical harmonic to use.
%> @param[in] N maximal degree of vector spherical harmonic to use.
%> @param[in] R Radius of the penetrable sphere.
%> @param[in] eps_r relative permittivity inside the sphere.
%> @param[out] F the far field operator of dimension [ (N+1)^2 - 1 ]^2
function F = ffoForPenetrableSphere(N, wavenum, R, eps_r)

Q = (N+1).^2 - 1;

F = zeros(2*Q,2*Q);

for n=1:N
    for m = -n:n
        
        indexU = n.^2 + n + m;
        indexV = indexU + Q;
        
        alpha_nm = zeros(1,Q);
        beta_nm = zeros(1,Q);
        beta_nm(indexU) = 4*pi*(1i)^n;

        Einc = fields.EntireWaveField( wavenum, alpha_nm, beta_nm );        
        Escat = scattering.scatFromPenetrableSphere(Einc, R, eps_r);
        [eInfPnm, eInfQnm] = Escat.farFieldCoeffs;
        
        F(1:Q,indexU) = eInfPnm;
        F(Q+1:2*Q,indexU) = eInfQnm;
                
        alpha_nm(indexU) = -4*pi*(1i)^n;
        beta_nm(indexU) = 0;

        Einc = fields.EntireWaveField( wavenum, alpha_nm, beta_nm );        
        Escat = scattering.scatFromPenetrableSphere(Einc, R, eps_r);
        [eInfPnm, eInfQnm] = Escat.farFieldCoeffs;
        
        F(1:Q,indexV) = eInfPnm;
        F(Q+1:2*Q,indexV) = eInfQnm;
        
    end
end


end

