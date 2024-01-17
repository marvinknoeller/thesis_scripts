%> @file L2tangFieldTest.m
%> @brief Tests for the L2tangFieldTest class
% ======================================================================
%> @brief Tests the functionality of the L2tangFieldTest class.
% ======================================================================
classdef L2tangFieldTest < matlab.unittest.TestCase
   
    methods (Test)
        
        %> @brief Test creation and evaluation from a RadiatingWaveField.
        %>
        %> A Hertz dipole is approximated by a RadiatingWaveField and the
        %> far field pattern of this function is represented as an L2tangField.
        %> The test then compares function evaluation of the far field of
        %> Hertz dipole with evaluations of the approximation.
        function testHertzDipoleApprox(testcase)
            
            N = 8;
            wavenum = sqrt(11);
            
            hd = fields.HertzDipole(wavenum, [0; 0.1; -0.1], ones(3,1)/sqrt(3));
            
            z =  [0;0;0];
            R = 2;
            
            rwf = fields.RadiatingWaveField.bestApprox(hd, wavenum, N, z, R); 
            
            [a, b] = rwf.farFieldCoeffs;
            farField = fields.L2tangField( a, b );
             
            % obtain the quadrature points and weights
            [Theta, Phi, W] = quadrature.sphereGaussLegTrap( 100, 100 );
            
            [F1,F2,F3] = farField.eval( Theta, Phi);
            [H1,H2,H3] = hd.farFieldPattern( Theta, Phi);
                        
            l2err = sqrt( abs(F1 - H1).^2 + abs(F2 - H2).^2 + abs(F3 - H3).^2 ) * W;
            fprintf('Error: %e\n\n',l2err);
            
            testcase.verifyLessThan(l2err,1e-5);
            
        end
        
        %> @brief Test the best approximation of a tangential vector field 
        %> which not a vector spherical harmonic by an L2tangField.
        %>
        %> The surface gradient of f(x) = exp( x1 + x2 + x3^2 ) is the function
        %> to be approximated.
        function testBestApproximation(testcase)
            
            N = 15;
            
            % define the surface gradient of f via a nested function and a wrapper field
            function [F1, F2, F3] = tangField(theta,phi)
                x1 = cos(phi).*sin(theta);
                x2 = sin(phi).*sin(theta);
                x3 = cos(theta);
                f = exp( x1 + x2 + x3.^2 );
                F1 = f .* (1 - x1.^2 - x1.*x2 - 2*x1.*x3.^2);
                F2 = f .* (1 - x1.*x2 - x2.^2 - 2*x2.*x3.^2);
                F3 = f .* (2*x3 - x1.*x3 - x2.*x3 - 2*x3.^3);
            end
            tf = fields.WrapperFieldOnSphere(@tangField);

            % Compute best approximation in L^2 norm
            l2tf = fields.L2tangField.bestApprox(tf,N);
            
            % obtain the quadrature points and weights
            [Theta, Phi, W] = quadrature.sphereGaussLegTrap( 100, 100 );
            [F1,F2,F3] = tf.eval(Theta,Phi);
            [G1,G2,G3] = l2tf.eval(Theta,Phi);
                        
            l2err = sqrt( abs(F1 - G1).^2 + abs(F2 - G2).^2 + abs(F3 - G3).^2 ) * W;
            fprintf('Error: %e\n\n',l2err);
            
            testcase.verifyLessThan(l2err,1e-5);
            
        end
    end
    
end