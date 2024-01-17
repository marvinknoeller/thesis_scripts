%> @file ScatteringTest.m
%> @brief Tests for the RadiatingWaveField class
% ======================================================================
%> @brief Tests the functionality of the functions in the package scattering.
% ======================================================================
classdef ScatteringTest < matlab.unittest.TestCase
   
    methods (Test)
        
        %> @brief Compute the scattered field for a perfectly conducting sphere and 
        %> check that the tangential components of the boundary values match.
        %>
        %> 
        function testScatFromPCSphere(testcase)
        
            % compute random wave number, center of asympotics and
            % incident wave field coefficients
            rng('shuffle');
            
            k = 5 * rand;
            y = rand(3,1);
            R = 1.5 + rand;
            
            N = 8;
            alpha_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            beta_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            
            Einc = fields.EntireWaveField( k, alpha_n, beta_n, y);
            
            Escat = scattering.scatFromPCSphere( Einc, R );                       
            
            % Evaluate both fields in 100 random points on the boundary of
            % the sphere and compare tangential components.       
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);      
            
            nu1 = cos(phi) .* sin(theta);
            nu2 = sin(phi) .* sin(theta);
            nu3 = cos(theta);
            
            X1 = y(1) + R .* nu1; 
            X2 = y(2) + R .* nu2; 
            X3 = y(3) + R .* nu3;
 
            [Ei1, Ei2, Ei3] = Einc.eval(X1,X2,X3);
            [Es1, Es2, Es3] = Escat.eval(X1,X2,X3);
            
            tangTot1 = nu2 .* (Ei3 + Es3) - nu3 .* (Ei2 + Es2);
            tangTot2 = nu3 .* (Ei1 + Es1) - nu1 .* (Ei3 + Es3);
            tangTot3 = nu1 .* (Ei2 + Es2) - nu2 .* (Ei1 + Es1);
            
            err = sqrt( abs(tangTot1).^2 + abs(tangTot2).^2 + abs(tangTot3).^2 );
            fprintf('Maximal error: %e\n', max(err) );
            
            testcase.verifyLessThan(err,1e-5);
            
        end
        
        %> @brief Compute the scattered field for a penetrable sphere and 
        %> check that the tangential components of the boundary values match.
        function testScatFromPenetrableSphere(testcase)
        
            % compute random wave number, relative permittivity, center 
            % of asympotics and incident wave field coefficients
            rng('shuffle');
            
            k = 5 * rand;
            eps_r = 2*rand -1 + 1i*rand;
            
            y = rand(3,1);
            R = 1.5 + rand;
            
            N = 3;
            alpha_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            beta_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            
            Einc = fields.EntireWaveField( k, alpha_n, beta_n, y);            
            [Escat, Etrans] = scattering.scatFromPenetrableSphere(Einc, R, eps_r);                       
            
            % Evaluate all three fields in 100 random points on the boundary of
            % the sphere and compute jump of tangential component.       
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1);      
            
            nu1 = cos(phi) .* sin(theta);
            nu2 = sin(phi) .* sin(theta);
            nu3 = cos(theta);
            
            X1 = y(1) + R .* nu1; 
            X2 = y(2) + R .* nu2; 
            X3 = y(3) + R .* nu3;
 
            [Ei1, Ei2, Ei3] = Einc.eval(X1,X2,X3);
            [Es1, Es2, Es3] = Escat.eval(X1,X2,X3);
            [Et1, Et2, Et3] = Etrans.eval(X1,X2,X3);
            
            tangTot1 = nu2 .* (Ei3 + Es3 - Et3) - nu3 .* (Ei2 + Es2 - Et2);
            tangTot2 = nu3 .* (Ei1 + Es1 - Et1) - nu1 .* (Ei3 + Es3 - Et3);
            tangTot3 = nu1 .* (Ei2 + Es2 - Et2) - nu2 .* (Ei1 + Es1 - Et1);
            
            err = sqrt( abs(tangTot1).^2 + abs(tangTot2).^2 + abs(tangTot3).^2 );
            fprintf('Maximal error in electric bc: %e\n', max(err) );
            
            testcase.verifyLessThan(err,1e-5);
            
            Hinc = Einc.getHFromE(1);
            Hscat = Escat.getHFromE(1);
            Htrans = Etrans.getHFromE(1/sqrt(eps_r));
 
            [Hi1, Hi2, Hi3] = Hinc.eval(X1,X2,X3);
            [Hs1, Hs2, Hs3] = Hscat.eval(X1,X2,X3);
            [Ht1, Ht2, Ht3] = Htrans.eval(X1,X2,X3);
            
            tangTot1 = nu2 .* (Hi3 + Hs3 - Ht3) - nu3 .* (Hi2 + Hs2 - Ht2);
            tangTot2 = nu3 .* (Hi1 + Hs1 - Ht1) - nu1 .* (Hi3 + Hs3 - Ht3);
            tangTot3 = nu1 .* (Hi2 + Hs2 - Ht2) - nu2 .* (Hi1 + Hs1 - Ht1);
            
            err = sqrt( abs(tangTot1).^2 + abs(tangTot2).^2 + abs(tangTot3).^2 );
            fprintf('Maximal error in magnetic bc: %e\n', max(err) );
            
            
        end
        
    end
    
end
