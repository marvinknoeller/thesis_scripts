%> @file EntireWaveFieldTest.m
%> @brief Tests for the EntireWaveField class
% ======================================================================
%> @brief Tests the functionality of the EntireWaveField class.
% ======================================================================
classdef EntireWaveFieldTest < matlab.unittest.TestCase
   
    methods (Test)
        
        
        %> @brief Compute and check accuracy of the evaluation of a plane wave.
        %>
        %> A plane wave is defined and its expansion coefficients as an
        %> EntireWaveFunction are computed. The result is checked for
        %> accuracy in 100 random points in a ball.
        function testEvaluation(testcase)
        
            % compute random wave number, direction of incidence and polarization
            rng('shuffle');
            
            z = rand(3,1);
            % R = 3;
            
            theta_d = pi * rand;
            phi_d = -pi + 2*pi*rand;
            d = [ cos(phi_d) * sin(theta_d); sin(phi_d) * sin(theta_d); cos(theta_d) ];
            
            q = -1 + 2*rand(3,1);
            p = cross(d,q);
            
            k = 5 * rand;
            % k = sqrt(11);
             
            %N = 25;
            N = 25;
            
            a_nm = zeros(1,(N+1)^2 - 1);
            b_nm = zeros(size(a_nm));

            [ d_phi, d_th, ~ ] = cart2sph( d(1), d(2), d(3) );
            d_th = pi/2 - d_th;

            for n=1:N

                m_range = n^2:((n+1)^2 - 1);

                [U_d, V_d] = specialfunc.vectorSphericalHarmonics( d_th, d_phi, n );
                p_times_U_d = p(1) * U_d(:,:,1)' + p(2) * U_d(:,:,2)' + p(3) * U_d(:,:,3)';
                p_times_V_d = p(1) * V_d(:,:,1)' + p(2) * V_d(:,:,2)' + p(3) * V_d(:,:,3)';
                
%                 % scaling factors
%                 j_n_over_kR = sqrt(pi / (2*k) / R ) * besselj(n+1/2,k*R) / (k*R);
%                 j_n_min_1_R = sqrt(pi / (2*k) / R ) .* besselj(n-1/2,k*R);
%                 scale_M_nm = k * R * j_n_over_kR;
%                 scale_curl_M_nm =  j_n_min_1_R - n * j_n_over_kR;
% 
%                 a_nm(m_range) = -4*pi * (1i)^n * p_times_V_d * exp( 1i * k * dot(z,d) ) * scale_M_nm;
%                 b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d * exp( 1i * k * dot(z,d) ) * scale_curl_M_nm;
% 
                a_nm(m_range) = -4*pi * (1i)^n * p_times_V_d * exp( 1i * k * dot(z,d) );
                b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d * exp( 1i * k * dot(z,d) );

            end
            
            pw  = fields.PlaneWave( k, d, p );
            ewf = fields.EntireWaveField(k, a_nm, b_nm, z);
            
            % Evaluate both fields in 100 random points in ball of radius 1.5 and compare results            
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);
            r = 1.5 * rand(100,1);            
            
            X1 = z(1) + r .* cos(phi) .* sin(theta); 
            X2 = z(2) + r .* sin(phi) .* sin(theta); 
            X3 = z(3) + r .* cos(theta);
            
            [PW1, PW2, PW3] = pw.eval(X1,X2,X3);
            [EWF1, EWF2, EWF3] = ewf.eval(X1,X2,X3);
            
            err = sqrt( abs(PW1 - EWF1).^2 + abs(PW2 - EWF2).^2 + abs(PW3 - EWF3).^2 );
            fprintf('Maximal error: %e\n', max(err) );
            
            testcase.verifyLessThan(err,1e-5);
            
%             % Evaluate both fields on the ball of radius 1.5 and compare
%             % results
%             
%             sphere = Sphere(z,1.5);
%             sphere.predefQuadPoints(50,50);
%             
%             [X1,X2,X3] = sphere.getCartesianQuadratureCoord;            
%             [PW1, PW2, PW3] = pw.eval(X1,X2,X3);
%             [EWF1, EWF2, EWF3] = ewf.evalOnSurf(sphere);
%             
%             err = sqrt( abs(PW1 - EWF1).^2 + abs(PW2 - EWF2).^2 + abs(PW3 - EWF3).^2 );
%             fprintf('Maximal error: %e\n', max(err) );
%             
%             testcase.verifyLessThan(err,1e-5);
            
        end
        
        %> @brief Compute and check accuracy of the approximation of a plane wave.
        %>
        %> A plane wave is defined and its best approximation by an
        %> EntireWaveFunction is computed. The result is checked for
        %> accuracy in 100 random points in a ball.
        function testBestApproximation(testcase)
            
            % Arguments for expansion
            z = [1; -0.2; 0.3];
            R = 2;
            N = 21;
        
            % compute random wave number, direction of incidence and polarization
            rng('shuffle');
            
            theta_d = pi * rand;
            phi_d = -pi + 2*pi*rand;
            d = [ cos(phi_d) * sin(theta_d); sin(phi_d) * sin(theta_d); cos(theta_d) ];
            
            q = -1 + 2*rand(3,1);
            p = cross(d,q);
            
            k = 5 * rand;

%             d = [ 1; 0; 0];
%             p = [ 0; 1; 0];
%             k = 4*sqrt(2);

            [ phi_d, theta_d, ~ ] = cart2sph( d(1), d(2), d(3) );
            theta_d = pi/2 - theta_d;
            
            % compute the exact coefficients of this plane wave.      
            a_nm = zeros(1,(N+1)^2 - 1);
            b_nm = zeros(size(a_nm));
            
            for n=1:N

                m_range = n^2:((n+1)^2 - 1);

                [U_d, V_d] = specialfunc.vectorSphericalHarmonics( theta_d, phi_d, n );
                p_times_U_d = p(1) * U_d(:,:,1)' + p(2) * U_d(:,:,2)' + p(3) * U_d(:,:,3)';
                p_times_V_d = p(1) * V_d(:,:,1)' + p(2) * V_d(:,:,2)' + p(3) * V_d(:,:,3)';
                
%                 % scaling factors
%                 j_n_over_kR = sqrt(pi / (2*k) / R ) * besselj(n+1/2,k*R) / (k*R);
%                 j_n_min_1_R = sqrt(pi / (2*k) / R ) .* besselj(n-1/2,k*R);
%                 scale_M_nm = k * R * j_n_over_kR;
%                 scale_curl_M_nm =  j_n_min_1_R - n * j_n_over_kR;
% 
%                 a_nm(m_range) = -4*pi * (1i)^n * p_times_V_d * exp( 1i * k * dot(z,d) ) * scale_M_nm;
%                 b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d * exp( 1i * k * dot(z,d) ) * scale_curl_M_nm;

                a_nm(m_range) = -4*pi * (1i)^n * p_times_V_d * exp( 1i * k * dot(z,d) );
                b_nm(m_range) =  4*pi * (1i)^n * p_times_U_d * exp( 1i * k * dot(z,d) );

            end
            
            % define plane wave
            pw = fields.PlaneWave( k, d, p );
            
            % compute EntireWaveField
            ewf = fields.EntireWaveField.bestApprox(pw, k, N, z, R);
            
%             A = [ a_nm; ewf.alpha; b_nm; ewf.beta ];
            
            % Evaluate both fields in 100 random points in ball of radius 1.5 and compare results            
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);
            r = 1.5 * rand(100,1);            
            
            X1 = z(1) + r .* cos(phi) .* sin(theta); 
            X2 = z(2) + r .* sin(phi) .* sin(theta); 
            X3 = z(3) + r .* cos(theta);

%             x = linspace(-1,1,101);
%             [x1,x2] = meshgrid(x,x);
%             x1 = z(1) + x1; x2 = z(2) + x2;
%             X1 = x1(:); X2 = x2(:); X3 = z(3) + zeros(size(X1));
            
            [PW1, PW2, PW3] = pw.eval(X1,X2,X3);
            [EWF1, EWF2, EWF3] = ewf.eval(X1,X2,X3);
            
            err = sqrt( abs(PW1 - EWF1).^2 + abs(PW2 - EWF2).^2 + abs(PW3 - EWF3).^2 );
            fprintf('Maximal error: %e\n', max(err) );
                       
%             pw2 = reshape(PW2,size(x1));
%             ewf2 = reshape(EWF2,size(x1));
%             
%             figure(1), surf(x1,x2,real(ewf2)), shading interp, view(2), colorbar
%             figure(2), surf(x1,x2,real(pw2)), shading interp, view(2), colorbar
            
            testcase.verifyLessThan(err,1e-5);
            
        end
        
    end
end