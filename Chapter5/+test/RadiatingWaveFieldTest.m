%> @file RadiatingWaveFieldTest.m
%> @brief Tests for the RadiatingWaveField class
% ======================================================================
%> @brief Tests the functionality of the RadiatingWaveField class.
% ======================================================================
classdef RadiatingWaveFieldTest < matlab.unittest.TestCase
   
    methods (Test)
        
        
        
%         %> A dipole is defined and its expansion coefficients as an
%         %> RadiatingWaveFunction are computed. The result is checked for
%         %> accuracy in 100 random points outside a ball.
%         function testEvaluation(testcase)
%         
%             % compute random wave number and source point
%             rng('shuffle');
%             
%             y = -1 + 2*rand(3,1);
%             k = 5 * rand;
%              
%             N = 25;
%             
%             % set center of expansion and minimal distance for x
%             z = [0; 0; 0];
%             [ yPhi, yTheta, yDist ] = cart2sph( y(1)-z(1), y(2)-z(2), y(3)-z(3) );
%             yTheta = pi/2 - yTheta;
%             Rmin = yDist + 1;
%             
%             a_nm = zeros(1,(N+1)^2 - 1);
%             b_nm = zeros(size(a_nm));
% 
%             for n=1:N
% 
%                 m_range = n^2:((n+1)^2 - 1);
%                 
%                 j_n  = sqrt(pi / (2*k) ./ yDist ) .* besselj(n+1/2,k*yDist);
%                 yFactor = 1i * k * sqrt(n*(n+1)) * j_n;
%                 Y_n = sphericalHarmonics( yTheta, yPhi, n);
% 
%                 a_nm(m_range) = yFactor * Y_n';
%                 %b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d;
% 
%             end
%             
%             pw = Dipole( k, y );
%             rwf = RadiatingWaveField(k, a_nm, b_nm, z, Rmin);
%             
%             % Evaluate both fields in 100 random points in ball of radius 1.5 and compare results            
%             theta = pi * rand(100,1);
%             phi = -pi + 2*pi*rand(100,1);
%             R = Rmin + 3 * rand(100,1);
%              
%             X1 = z(1) + R .* cos(phi) .* sin(theta); 
%             X2 = z(2) + R .* sin(phi) .* sin(theta); 
%             X3 = z(3) + R .* cos(theta);
%             
%             [PW1, PW2, PW3] = pw.eval(X1,X2,X3);
%             [RWF1, RWF2, RWF3] = rwf.eval(X1,X2,X3);
%             
%             err = sqrt( abs(PW1 - RWF1).^2 + abs(PW2 - RWF2).^2 + abs(PW3 - RWF3).^2 );
%             fprintf('Maximal error: %e\n', max(err) );
%             
%             testcase.verifyLessThan(err,1e-5);
%             
%             
%         end
        
        %> @brief Compute and check accuracy of the approximation of a
        %> Hertz dipole.
        %>
        %> A Hertz dipole is defined and its best 
        %> approximation by a RadiatingWaveFunction is computed. The result 
        %> is checked for accuracy in 100 random points outside the sphere
        %> on which the best approximation is computed, and again on a
        %> sphere with the same center but a larger radius.
        function testBestApproximation(testcase)
            
            % Arguments for expansion
            z = [0.0; 0.0; 0.5];
            R = 3;
            N = 19;
        
            % set random wave number, source point and polarization for the
            % Hertz dipole
            rng('shuffle');
            
            y = -1 + 2*rand(3,1);
            k = 5 * rand;
            p_angles = [0;-pi] + pi*[1;2] .* rand(2,1);
            p = [ cos(p_angles(2)) * sin(p_angles(1)); sin(p_angles(2))*sin(p_angles(1)); cos(p_angles(1)) ];
            
            % define the dipole           
            hd = fields.HertzDipole( k, y, p );
            
            % compute RadiatingWaveField
            rwf = fields.RadiatingWaveField.bestApprox(hd, k, N, z, R); 
            rwf.printNormInDegree();
            
            % Evaluate both fields in 100 random points in a shell of minimal 
            % radius 3 and maximal radius 6 and compare results            
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);
            r = R +  3 * rand(100,1);            
            
            X1 = z(1) + r .* cos(phi) .* sin(theta); 
            X2 = z(2) + r .* sin(phi) .* sin(theta); 
            X3 = z(3) + r .* cos(theta);
 
            [HD1, HD2, HD3] = hd.eval(X1,X2,X3);
            [RWF1, RWF2, RWF3] = rwf.eval(X1,X2,X3);
            
            err = sqrt( abs(HD1 - RWF1).^2 + abs(HD2 - RWF2).^2 + abs(HD3 - RWF3).^2 );
            fprintf('Maximal error: %e\n', max(err) );
            
            testcase.verifyLessThan(err,1e-5);
            
        end
        
        %> @brief Compute and check accuracy for evaluations of the far fields
        %> of first degree radiating wave fields.
        %>
        %> The six first order radiating wave fields (Order 1 spherical Hankel 
        %> function times first degree spherical vector harmonic) for random
        %> centers of expansion are approximated on a sphere around the
        %> origin, the far field coefficients of the approximation are
        %> computed and the far field is evaluated in 100 random points on
        %> the unit sphere. The results are compared to the exact values.
        function testFarFieldEvaluation(testcase)
                                
            N = 15;
            
            % compute random wave number and center of expansion
            rng('shuffle');
            
            y = -1 + 2*rand(3,1);
            k = 5 * rand;
            
%             y = [0;0;0];
%             k = 3*sqrt(2);

            % Define the wave fields
            rwf = cell(1,6);
            for j=1:3
                alpha = zeros(1,3);
                alpha(j) = 1;
                rwf{j} = fields.RadiatingWaveField( k, alpha, zeros(1,3), y );
                beta = zeros(1,3);
                beta(j) = 1;
                rwf{j+3} = fields.RadiatingWaveField( k, zeros(1,3), beta, y );
            end
            
            % Compute best approximations and their far fields
            R = 4;
            z = [0;0;0];
            bestApprox = cell(1,6);
            bestApproxFF = cell(1,6);
            
            for j=1:6
                bestApprox{j} = fields.RadiatingWaveField.bestApprox( rwf{j}, k, N, z, R);
                [a, b] = bestApprox{j}.farFieldCoeffs;
                bestApproxFF{j} = fields.L2tangField( a, b );
            end

%              t = linspace(-2,2,101);
%              [X1,X2] = meshgrid(R + 3 + t, t);
%              X3 = zeros(size(X1));
%  
%             [PW1, PW2, PW3] = bestApprox{5}.eval(X1,X2,X3);
%             [RWF1, RWF2, RWF3] = rwf{5}.eval(X1,X2,X3);
%              
%              figure(1), surf(X1,X2,real(PW3)), shading interp, view(2), axis equal, colorbar
%              figure(2), surf(X1,X2,real(RWF3)), shading interp, view(2), axis equal, colorbar
            
            farFieldError = zeros(100,6);
            farFieldError2 = zeros(100,6);
            
            % Evaluate far fields in 100 random points on unit sphere and compare results            
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);
            
            % Phase factor to compensate for center of expansion not in the
            % origin.
            factorCenterExp = exp( -1i * k * ( cos(phi) .* sin(theta) * y(1) + sin(phi) .* sin(theta) * y(2) + cos(theta) * y(3) ) );
                
%             % scaling factors
%             h_1_over_k  = sqrt(pi / (2*k) ) .* besselh(3/2,1,k) ./ k;
%             h_0_1 = sqrt(pi / (2*k) ) .* besselh(1/2,1,k);
%             scale_N_1m = k * h_1_over_k;
%             scale_curl_N_1m =  h_0_1 - h_1_over_k;
            
            for j=1:3
                [U, V] = specialfunc.vectorSphericalHarmonics( theta.', phi.', 1 );
                [W1, W2, W3] = bestApproxFF{j}.eval(theta,phi);
                farFieldError(:,j) = sqrt( ...
                    abs( W1 - 4*pi/k * factorCenterExp .* V(j,:,1).' ).^2 ...
                    + abs( W2 - 4*pi/k * factorCenterExp .* V(j,:,2).' ).^2 ...
                    + abs( W3 - 4*pi/k * factorCenterExp .* V(j,:,3).' ).^2 ...
                    );
                [Z1, Z2, Z3] = rwf{j}.evalFarFieldPattern(theta,phi);
                farFieldError2(:,j) = sqrt( ...
                    abs( Z1 - 4*pi/k * factorCenterExp .* V(j,:,1).' ).^2 ...
                    + abs( Z2 - 4*pi/k * factorCenterExp .* V(j,:,2).' ).^2 ...
                    + abs( Z3 - 4*pi/k * factorCenterExp .* V(j,:,3).' ).^2 ...
                    );
                
                [W1, W2, W3] = bestApproxFF{j+3}.eval(theta,phi);
                farFieldError(:,j+3) = sqrt( ...
                    abs( W1 + 4*pi/k * factorCenterExp .* U(j,:,1).' ).^2 ...
                    + abs( W2 + 4*pi/k * factorCenterExp .* U(j,:,2).' ).^2 ...
                    + abs( W3 + 4*pi/k * factorCenterExp .* U(j,:,3).' ).^2 ...
                    );
                [Z1, Z2, Z3] = rwf{j+3}.evalFarFieldPattern(theta,phi);
                farFieldError(:,j+3) = sqrt( ...
                    abs( Z1 + 4*pi/k * factorCenterExp .* U(j,:,1).' ).^2 ...
                    + abs( Z2 + 4*pi/k * factorCenterExp .* U(j,:,2).' ).^2 ...
                    + abs( Z3 + 4*pi/k * factorCenterExp .* U(j,:,3).' ).^2 ...
                    );
            end
            
            fprintf('Maximum Far Field Error: %g\n', max(farFieldError));
            fprintf('Maximum Far Field Error direct evaluation: %g\n', max(farFieldError2));
            
            testcase.verifyLessThan(max(max(farFieldError)),1e-5);
            
        end
    end
    
end
