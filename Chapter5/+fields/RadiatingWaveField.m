%> @file RadiatingWaveField.m
%> @brief Contains the RadiatingWaveField class.
% ======================================================================
%> @brief Class representing a radiating vector field defined outside a
%> given ball.
%>
%> Radiating wave functions are linear combinations of the functions
%>
%> @f[
%> N_n^m(x) = - h_n^{(1)}(k |x|) \, V_n^m(\hat{x}) 
%> \qquad \text{and} \qquad
%> \frac{1}{i k} \, \operatorname{curl} \, N_n^m(x)
%> @f]
%>
%> which are solutions to the Maxwell system. More precisely, we shift
%> these functions to have an arbitrary center of expansion z. Hence, an 
%> object of this class represents a function
%>
%> @f[
%> W(x) = \sum_{n=1}^N \sum_{m = -n}^n \left[ \alpha_n^m \, N_n^m (x-z) +
%> \beta_n^m \, \frac{1}{i k} \, \operatorname{curl} \, N_n^m(x-z) \right]
%> @f]
%>
%> Typically, the coefficients are computed by a best approximation on a
%> sphere of radius R around z. Then
%>
%> @f[
%>  \alpha_n^m = - \int_{S^2} \frac{ W(z + R \hat{x}) \cdot V_n^{-m}(\hat{x}) }{ h_n^{(1)}(kR) } \,
%>   ds(\hat{x}) \, , 
%>   \qquad 
%>   \beta_n^m = i \int_{S^2} \frac{ W(z + R \hat{x}) \cdot U_n^{-m}(\hat{x}) }{
%>     \, h_{n-1}^{(1)}(kR) - n/(kR) \, h_n^{(1)}(kR) } \, ds(\hat{x})
%> @f]
%
% ======================================================================
classdef RadiatingWaveField < fields.VectorField
    
    properties (SetAccess = protected)
        
        %> @brief The wave number
        wavenum = 1;
        
        %> @brief The degree of the highest vector spherical harmonic in
        %> the expansion.
        degree = 1;
        
        %> @brief The coordinates of the center of expansion, a 3
        %> coordinate column vector.
        z = [0; 0; 0];
        
        %> @brief The coefficients for the functions @f$ N_n^{m} @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like
        %>
        %> a_1^{-1} a_1^{0} a_1^{1} a_2^{-2} ... a_2^{2} a_3^{-3} ... a_degree^{degree}
        %>
        alpha = [];
        
        %> @brief The coefficients for the functions @f$ \frac{1}{i k} \, \operatorname{curl} \, N_n^m @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like the coefficients in @a alpha.
        beta = [];
        
    end
    
    
    methods
        
        %> @brief Construction
        %>
        %> Precise behaviour depends on the number of arguments:
        %> - at least 1 argument:  field around the origin
        %> - at least 4 arguments: field around z
        %>
        %> If at least 3 arguments are given, the coefficients are set.
        %> Otheriwse, all coeffients are initialized to 0.
        function rwf = RadiatingWaveField( wavenum, alpha, beta, z )
            
            
            if ( nargin >= 1 )
                rwf.wavenum = wavenum;
            end
            if ( nargin >= 3 )
                rwf.setCoeffs( alpha, beta );                
            end
            
            if ( nargin >= 4 )
                rwf.z = z;
            end
            
        end
        
        %> @brief Set center of expansion and minimal radius
        function setCenterRad(rwf, z)
            
            %> @todo Check admissability of arguments
            
            rwf.z = z;
            
        end
        
        %> @brief Set degree and expanson coefficients
        function setCoeffs(rwf, alpha, beta)
            
            %> @todo Check sizes
            
            rwf.alpha = alpha;
            rwf.beta = beta;
            rwf.degree = sqrt( length(alpha) + 1 ) - 1;            
            
        end
        
        %> @brief Implementation of VectorField.eval()
        function [V1, V2, V3] = eval( rwf, X1, X2, X3 )
            
            x1 = X1(:).' - rwf.z(1);
            x2 = X2(:).' - rwf.z(2);
            x3 = X3(:).' - rwf.z(3);

            V1 = zeros(1,length(x1));
            V2 = zeros(1,length(x1));
            V3 = zeros(1,length(x1));

            [ phi, theta, r ] = cart2sph( x1, x2, x3 );
            theta = pi/2 - theta;

            for n=1:rwf.degree

                % Compute Hankel functions
                h_n_over_kr  = sqrt(pi / (2*rwf.wavenum) ./ r ) .* besselh(n+1/2,1,rwf.wavenum*r) ./ ( rwf.wavenum*r );
                h_n_min_1_r = sqrt(pi / (2*rwf.wavenum) ./ r ) .* besselh(n-1/2,1,rwf.wavenum*r);
%                 h_n_over_kR  = sqrt(pi / (2*rwf.wavenum) ./ rwf.sphere.R ) .* besselh(n+1/2,1,rwf.wavenum*rwf.sphere.R) ./ ( rwf.wavenum*rwf.sphere.R );
%                 h_n_min_1_R = sqrt(pi / (2*rwf.wavenum) ./ rwf.sphere.R ) .* besselh(n-1/2,1,rwf.wavenum*rwf.sphere.R);
%                 
%                 % scaling factors
%                 scale_N_nm = rwf.wavenum * rwf.sphere.R * h_n_over_kR;
%                 scale_curl_N_nm =  h_n_min_1_R - n * h_n_over_kR;
                
                % despite the variable names, these are only the factors
                % that only depend on the radial variable
                curl_N_nm_Y = -1i * sqrt(n*(n+1))  * h_n_over_kr;
                curl_N_nm_U = -1i * ( h_n_min_1_r - n * h_n_over_kr );
                N_nm = -rwf.wavenum * r .* h_n_over_kr;

                Y_x = specialfunc.sphericalHarmonics( theta, phi, n);
                [U_x, V_x] = specialfunc.vectorSphericalHarmonics( theta, phi, n );
                
                V1 = V1 + curl_N_nm_Y .* ( rwf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * ( cos(phi) .* sin(theta) ) ) ) ) ...
                    + curl_N_nm_U .* ( rwf.beta(n^2:((n+1)^2-1)) * U_x(:,:,1) ) + N_nm .* ( rwf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,1) );
                V2 = V2 + curl_N_nm_Y .* ( rwf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * ( sin(phi) .* sin(theta) ) ) ) ) ...
                    + curl_N_nm_U .* ( rwf.beta(n^2:((n+1)^2-1)) * U_x(:,:,2) ) + N_nm .* ( rwf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,2) );
                V3 = V3 + curl_N_nm_Y .* ( rwf.beta(n^2:((n+1)^2-1)) * ( Y_x .* ( ones(2*n+1,1) * cos(theta) ) ) ) ...
                    + curl_N_nm_U .* ( rwf.beta(n^2:((n+1)^2-1)) * U_x(:,:,3) ) + N_nm .* ( rwf.alpha(n^2:((n+1)^2-1)) * V_x(:,:,3) );

            end
            
            V1 = reshape( V1, size(X1) );
            V2 = reshape( V2, size(X1) );
            V3 = reshape( V3, size(X1) );

        end
        
        
        %> @brief Coefficients of the far field pattern of the RadiatingWaveField
        %>
        %> F^\infty = \sum a_n^m U_n^m + b_n^m V_n^m
        %>
        %> This function computes the coefficients for the far field
        %> asympotics using z as the origin.
        function [a, b] = farFieldCoeffs(rwf)
            
            preFac_alpha = ones(1,(rwf.degree+1)^2-1);
            preFac_beta = ones(1,(rwf.degree+1)^2-1);
            
            for n=1:rwf.degree
                mRange = n^2:( (n+1)^2 - 1 );
                
%                 % scaling factors
%                 h_n_over_kR  = sqrt(pi / (2*rwf.wavenum) ./ rwf.sphere.R ) .* besselh(n+1/2,1,rwf.wavenum*rwf.sphere.R) ./ ( rwf.wavenum*rwf.sphere.R );
%                 h_n_min_1_R = sqrt(pi / (2*rwf.wavenum) ./ rwf.sphere.R ) .* besselh(n-1/2,1,rwf.wavenum*rwf.sphere.R);
%                 scale_N_nm = rwf.wavenum * rwf.sphere.R * h_n_over_kR;
%                 scale_curl_N_nm =  h_n_min_1_R - n * h_n_over_kR;
%                 
%                 preFac_alpha(mRange) = -4*pi/rwf.wavenum * (-1i)^(n+1) / scale_N_nm;
%                 preFac_beta(mRange) = 4*pi/rwf.wavenum * (-1i)^(n+1) / scale_curl_N_nm;

                preFac_alpha(mRange) = -4*pi/rwf.wavenum * (-1i)^(n+1);
                preFac_beta(mRange) = 4*pi/rwf.wavenum * (-1i)^(n+1);
            end
            
            a = preFac_beta .* rwf.beta;
            b = preFac_alpha .* rwf.alpha;
            
        end
        
        %> @brief Evaluation of the far field pattern of this wave field
        %>
        %> For the far field asympotics, the point y is used as the origin.
        %> If this argument is omitted, y = [0;0;0] is assumed.
        %>
        %> The far field pattern is evaluated for the directions given by 
        %> @a theta and @a phi in angular coordinates. This is done by
        %> obtaining the far field as an L2TangField if z were at the origin. 
        %> This pattern is evaluated and then multiplied by the appropriate 
        %> phase shift.
        function [V1, V2, V3] = evalFarFieldPattern( rwf, theta, phi, y )
            
            if (nargin < 4)
                y = [0;0;0];
            end
            
            [a, b] = rwf.farFieldCoeffs;
            l2t = fields.L2tangField( a, b );
            
            [V1, V2, V3] = l2t.eval( theta, phi );
            
            phaseFactor = exp( -1i * rwf.wavenum ...
                * ( cos(phi) .* sin(theta) * ( rwf.z(1) - y(1) ) ...
                  + sin(phi) .* sin(theta) * ( rwf.z(2) - y(2) ) ...
                  + cos(theta) * ( rwf.z(3) - y(3) ) ) );            
              
            V1 = phaseFactor .* V1;
            V2 = phaseFactor .* V2;
            V3 = phaseFactor .* V3;
            
        end % of function [V1, V2, V3] = evalFarFieldPattern( rwf, theta, phi, y )
        
                        
        %> @brief Obtain the far field pattern of this wave field as a
        %> fields.L2tangField
        %>
        %> The best approximation of the far field by a fields.L2tangField 
        %> of degree N is obtained with y used as the origin for the far field
        %> asympotics.        
        function l2t = obtainFarField( rwf, N, y)
            
            evalFunction = @(Theta, Phi) rwf.evalFarFieldPattern(Theta, Phi, y);            
            wf = fields.WrapperFieldOnSphere( evalFunction );
            
            l2t = fields.L2tangField.bestApprox( wf, N );
            
        end
        
        
        %> @brief Print norm in each degree
        function printNormInDegree(rwf)
            
            for n = 1:rwf.degree
                
                mRange = n^2:((n+1)^2-1);
                fprintf('  %10.4e', sqrt( sum( abs(rwf.alpha(mRange)).^2 + abs(rwf.beta(mRange)).^2 ) ) );
            end
            
            fprintf('\n');
            
        end
        
        %> @brief Obtain the represenation of the magnetic field from the
        %> electric field.
        %>
        %> The RadiatingWaveField is intepreted to be an electric field. Using
        %> the Maxwell system, the corresponding magnetic field is computed.
        %> This requires the impedance @f$ Z = \sqrt{ \mu / \varepsilon } @f$
        %> as an additional parameter
        function H = getHFromE(rwf, Z)
            
            a = -rwf.beta ./ Z;
            b = rwf.alpha ./ Z;
            
            H = fields.RadiatingWaveField( rwf.wavenum, a, b, rwf.z );
        end
        
        %> @brief Obtain the represenation of the electric field from the
        %> magnetic field.
        %>
        %> The RadiatingWaveField is intepreted to be an magnetic field. Using
        %> the Maxwell system, the corresponding electric field is computed.
        %> This requires the impedance @f$ Z = \sqrt{ \mu / \varepsilon } @f$
        %> as an additional parameter
        function E = getEFromH(rwf, Z)
            
            a = Z .* rwf.beta;
            b = -Z .* rwf.alpha;
            
            E = fields.RadiatingWaveField( rwf.wavenum, a, b, rwf.z );
        end
                       
    end % of methods
    
    
    methods (Static)
        
        %> @brief Compute the best approximation of a given VectorField by a RadiatingWaveField
        %>
        %> The VectorField @a vf is approximated by a RadiatingWaveField
        %> outside the ball of radius R centered at z. This done by computing 
        %> the best approximation to vf in the space spanned by the basis functions of 
        %> degree less than or equal to @a N on the sphere of radius R centered 
        %> at z  in the L^2-norm (i.e. compute an orthogonal projection)
        %>
        %> The function performs integration on the unit sphere with a
        %> tensor product quadrature rule, with a Gauss-Legendre rule in the 
        %> polar angle and a composite trapezoidal rule in the azimuthal angle. 
        %> 2*Nazim points in the
        %> azimuthal and Npolar points in the polar direction are used. If these
        %> arguments are not specified, a default value of 100 points in polar
        %> and 200 points in azimuthal direction is used.
        function rwf = bestApprox(vf, wavenum, N, z, R, Npolar, Nazim)
                    
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
                h_n_over_kR = sqrt(pi / (2*wavenum) / R ) * besselh(n+1/2,1,wavenum*R) / (wavenum*R);
                h_n_min_1_R = sqrt(pi / (2*wavenum) ./ R ) .* besselh(n-1/2,1,wavenum*R);
                scale_N_nm = wavenum * R * h_n_over_kR;
                scale_curl_N_nm =  h_n_min_1_R - n * h_n_over_kR;
                
                Wmat = W * ones(1,2*n+1);
                
                % evaluate the vsh on the quadrature points
                [U, V] = specialfunc.vectorSphericalHarmonics( Theta, Phi, n);
                
                % scalar producs in L^2_tangential
                a = F1 * ( V(:,:,1)' .* Wmat ) + F2 * ( V(:,:,2)' .* Wmat ) + F3 * ( V(:,:,3)' .* Wmat );
                b = F1 * ( U(:,:,1)' .* Wmat ) + F2 * ( U(:,:,2)' .* Wmat ) + F3 * ( U(:,:,3)' .* Wmat );

                % set the coefficients
                index = n^2:( (n+1)^2 - 1 );
                alpha(index) = -a / scale_N_nm;
                beta(index) = 1i * b / scale_curl_N_nm;
                
            end
            
            % define the RadiatingWaveField
            rwf = fields.RadiatingWaveField;
            rwf.wavenum = wavenum;
            rwf.z = z;
            rwf.setCoeffs( alpha, beta );
            
        end
        
        
    end % of methods (Static)
    
end
