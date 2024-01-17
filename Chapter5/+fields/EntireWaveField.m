%> @file EntireWaveField.m
%> @brief Contains the EntireWaveField class.
% ======================================================================
%> @brief Class representing a vector field that is an entire function
%> in a given ball.
%>
%> Entire wave functions are linear combinations of the functions
%>
%> @f[
%> M_n^m(x) = - j_n(k |x|) \, V_n^m(\hat{x}) 
%> \qquad \text{and} \qquad
%> \frac{1}{i k} \, \operatorname{curl} \, M_n^m(x)
%> @f]
%>
%> which are solutions to the Maxwell system. More precisely, we shift
%> these functions to have an arbitrary center of expansion z. Hence, an object of this
%> class represents a function
%>
%> @f[
%> W(x) = \sum_{n=1}^N \sum_{m = -n}^n \left[ \alpha_n^m \, M_n^m (x-z) +
%> \beta_n^m \, \frac{1}{i k} \, \operatorname{curl} \, M_n^m(x-z) \right]
%> @f]
%>
%> Typically, the coefficients are computed by a best approximation on a
%> sphere of radius R around z. Then
%>
%> @f[
%>  \alpha_n^m = - \int_{S^2} \frac{ W(z + R \hat{x}) \cdot V_n^{-m}(\hat{x}) }{ j_n(kR) } \,
%>   ds(\hat{x}) \, , 
%>   \qquad 
%>   \beta_n^m = i \int_{S^2} \frac{ W(z + R \hat{x}) \cdot U_n^{-m}(\hat{x}) }{
%>     \, j_{n-1}(kR) - n/(kR) \, j_n(kR) } \, ds(\hat{x})
%> @f]
%
% ======================================================================
classdef EntireWaveField < fields.VectorField
    % classdef EntireWaveField < VectorField
    
    properties (SetAccess = protected)
        
        %> @brief The wave number
        wavenum = 1;
        
        %> @brief The degree of the highest vector spherical harmonic in
        %> the expansion.
        degree = 1;
        
        %> @brief The coordinates of the center of expansion, a 3
        %> coordinate column vector.
        z = [0; 0; 0]
        
        %> @brief The coefficients for the functions @f$ M_n^{m} @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like
        %>
        %> a_1^{-1} a_1^{0} a_1^{1} a_2^{-2} ... a_2^{2} a_3^{-3} ... a_degree^{degree}
        %>
        alpha = [];
        
        %> @brief The coefficients for the functions @f$ \frac{1}{i k} \, \operatorname{curl} \, M_n^m @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like the coefficients in @a alpha.
        beta = [];
        
    end
    
    methods
        
        %> @brief Construction
        %>
        %> Precise behaviour depends on the number of arguments:
        %> - 3 arguments: field centered at origin
        %> - 4 arguments: field centered at z
        function ewf = EntireWaveField( wavenum, alpha, beta, z )
            
            if ( nargin >= 1 )
                ewf.wavenum = wavenum;
            end
            if ( nargin >= 3 )
                ewf.setCoeffs( alpha, beta );                
            end
            
            if ( nargin >= 4 )
                ewf.z = z;
            end
            
        end
        
        %> @brief Set center of expansion 
        function setCenter(ewf, z)
            
            %> @todo Check admissability of arguments
            
            ewf.z = z;
            
        end
        
        %> @brief Set degree and expanson coefficients
        function setCoeffs(ewf, alpha, beta)
            
            %> @todo Check sizes
            
            ewf.alpha = alpha;
            ewf.beta = beta;
            ewf.degree = sqrt( length(alpha) + 1 ) - 1;            
            
        end
        
        %> @brief Implementation of VectorField.eval()
        %>
        %> This function computes everything from scratch. No stored
        %> information is used.
        function [V1, V2, V3] = eval( ewf, X1, X2, X3 )
            
            x1 = X1(:).' - ewf.z(1);
            x2 = X2(:).' - ewf.z(2);
            x3 = X3(:).' - ewf.z(3);

            V1 = zeros(1,length(x1));
            V2 = zeros(1,length(x1));
            V3 = zeros(1,length(x1));

            [ phi, theta, r ] = cart2sph( x1, x2, x3 );
            theta = pi/2 - theta;

            for n=1:ewf.degree

                % Intitialize variables
                j_n_over_kr = zeros(size(r));
                j_n_min_1 = zeros(size(r));

                % for small arguments, we have to use asymptotic forms for the
                % Bessel functions
                r_small = ( r < 1e-5 );
                j_n_over_kr(r_small) = sqrt(pi / 16 ) /  gamma(n+3/2) * ( ewf.wavenum * r(r_small)/2 ).^(n-1);
                j_n_over_kr(~r_small)  = sqrt(pi / (2*ewf.wavenum) ./ r(~r_small) ) .* besselj(n+1/2,ewf.wavenum*r(~r_small)) ./ (ewf.wavenum*r(~r_small));
                j_n_min_1(r_small) = sqrt(pi / 4 ) / gamma(n+1/2) .* ( ewf.wavenum*r(r_small)/2 ).^(n-1);
                j_n_min_1(~r_small) = sqrt(pi / (2*ewf.wavenum) ./ r(~r_small) ) .* besselj(n-1/2,ewf.wavenum*r(~r_small));
                
                % despite the variable names, these are in fact the factors
                % that only depend on the radial variable
                curl_M_nm_Y = -1i * sqrt(n*(n+1)) * j_n_over_kr;
                curl_M_nm_U = -1i * ( j_n_min_1 - n * j_n_over_kr );
                M_nm = -ewf.wavenum * r .* j_n_over_kr;

                Y_x = specialfunc.sphericalHarmonics( theta, phi, n);
                [U_x, V_x] = specialfunc.vectorSphericalHarmonics( theta, phi, n );
                
                V1 = V1 + curl_M_nm_Y .* ( ewf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * ( cos(phi) .* sin(theta) ) ) ) ) ...
                    + curl_M_nm_U .* ( ewf.beta(n^2:((n+1)^2-1)) * U_x(:,:,1) ) + M_nm .* ( ewf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,1) );
                V2 = V2 + curl_M_nm_Y .* ( ewf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * ( sin(phi) .* sin(theta) ) ) ) ) ...
                    + curl_M_nm_U .* ( ewf.beta(n^2:((n+1)^2-1)) * U_x(:,:,2) ) + M_nm .* ( ewf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,2) );
                V3 = V3 + curl_M_nm_Y .* ( ewf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * cos(theta) ) ) ) ...
                    + curl_M_nm_U .* ( ewf.beta(n^2:((n+1)^2-1)) * U_x(:,:,3) ) + M_nm .* ( ewf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,3) );

            end
            
            % Reshape output to size of the input
            V1 = reshape( V1, size(X1) );
            V2 = reshape( V2, size(X1) );
            V3 = reshape( V3, size(X1) );

        end
        
        
        %> @brief Print norm in each degree
        %>
        %> Function for debugging purposes that prints the norm of
        %> projections onto functions of just one degree n.
        function printNormInDegree(ewf)
            
            for n = 1:ewf.degree
                
                mRange = n^2:((n+1)^2-1);
                fprintf('  %10.4e', sqrt( sum( abs(ewf.alpha(mRange)).^2 + abs(ewf.beta(mRange)).^2 ) ) );
            end
            
            fprintf('\n');
            
        end
        
        %> @brief Obtain the represenation of the magnetic field from the
        %> electric field.
        %>
        %> The EntireWaveField is intepreted to be an electric field. Using
        %> the Maxwell system, the corresponding magnetic field is computed.
        %> This requires the impedance @f$ Z = \sqrt{ \mu / \varepsilon } @f$
        %> as an additional parameter
        function H = getHFromE(ewf, Z)
            
            a = -ewf.beta ./ Z;
            b = ewf.alpha ./ Z;
            
            H = fields.EntireWaveField( ewf.wavenum, a, b, ewf.z );
        end
        
        %> @brief Obtain the represenation of the electric field from the
        %> magnetic field.
        %>
        %> The EntireWaveField is intepreted to be an magnetic field. Using
        %> the Maxwell system, the corresponding electric field is computed.
        %> This requires the impedance @f$ Z = \sqrt{ \mu / \varepsilon } @f$
        %> as an additional parameter
        function E = getEFromH(ewf, Z)
            
            a = Z .* ewf.beta;
            b = -Z .* ewf.alpha;
            
            E = fields.EntireWaveField( ewf.wavenum, a, b, ewf.z );
        end
        
    end
    
    methods (Static)
        
        %> @brief Compute the best approximation of a given VectorField by
        %> an EntireWaveField.
        %>
        %> The VectorField @a vf is approximated by an EntireWaveField
        %> within a ball of radius R centred at z. This done by computing the 
        %> best approximation  to vf in the space spanned by the basis functions 
        %> of degree less than or equal to @a N on the boundary of the ball 
        %> in the L^2-norm (i.e. compute an orthogonal projection)
        %>
        %> The function performs integration on the unit sphere with a
        %> tensor product quadrature rule, with a Gauss-Legendre rule in the 
        %> polar angle and a composite trapezoidal rule in the azimuthal angle. 
        %> 2*Nazim points in the
        %> azimuthal and Npolar points in the polar direction are used. If these
        %> arguments are not specified, a default value of 100 points in polar
        %> and 200 points in azimuthal direction is used.
        function ewf = bestApprox(vf, wavenum, N, z, R, Npolar, Nazim)
                    
            % obtain the quadrature points and weights
            if ( nargin >= 6 )
                [Theta, Phi, W] = quadrature.sphereGaussLegTrap( Npolar, Nazim );
            else
                [Theta, Phi, W] = quadrature.sphereGaussLegTrap( 100, 100 );
            end         
            
            X1 = z(1) + R * cos(Phi) .* sin(Theta);
            X2 = z(2) + R * sin(Phi) .* sin(Theta);
            X3 = z(3) + R * cos(Theta);
            
            % evaluate the VectorField in the quadrature points
            [ F1, F2, F3 ] = vf.eval( X1, X2, X3 );
            
            % intialize coefficients
            alpha = zeros(1,(N+1)^2 - 1);
            beta = zeros( size(alpha) );
            
            % for each vector spherical harmonic up to degree N
            for n=1:N
                
                % scaling of coefficients on sphere of radius R
                j_n_over_kR = sqrt(pi / (2*wavenum) / R ) * besselj(n+1/2,wavenum*R) / (wavenum*R);
                j_n_min_1_R = sqrt(pi / (2*wavenum) ./ R ) .* besselj(n-1/2,wavenum*R);
                scale_M_nm = wavenum * R * j_n_over_kR;
                scale_curl_M_nm =  j_n_min_1_R - n * j_n_over_kR;
                
                Wmat = W * ones(1,2*n+1);
                
                % evaluate the vsh on the quadrature points
                [U, V] = specialfunc.vectorSphericalHarmonics( Theta, Phi, n);
                
                % scalar producs in L^2_tangential
                a = F1 * ( V(:,:,1)' .* Wmat ) + F2 * ( V(:,:,2)' .* Wmat ) + F3 * ( V(:,:,3)' .* Wmat );
                b = F1 * ( U(:,:,1)' .* Wmat ) + F2 * ( U(:,:,2)' .* Wmat ) + F3 * ( U(:,:,3)' .* Wmat );

                % set the coefficients
                index = n^2:( (n+1)^2 - 1 );
                alpha(index) = -a / scale_M_nm;
                beta(index) = 1i * b / scale_curl_M_nm;
                
            end
            
            % define the EntireWaveField
            ewf = fields.EntireWaveField(wavenum);
            ewf.setCenter(z);
            ewf.setCoeffs( alpha, beta );
            
        end
        
        
    end % of methods (Static)
    
end

