%> @file WrapperFieldOnSphereTest.m
%> @brief Tests for the WrapperFieldOnSphereTest class
% ======================================================================
%> @brief Tests the functionality of the WrapperFieldOnSphere class.
% ======================================================================
classdef WrapperFieldOnSphereTest < matlab.unittest.TestCase
   
    methods (Test)
        
        function testEvalCorrect(testcase)
        
            % compute random wave number, center of asympotics and
            % radiating wave field coefficients
            rng('shuffle');
            
            k = 5 * rand;
            y = rand(3,1);
            
            N = 8;
            alpha_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            beta_n = rand(1,(N+1)^2 - 1) + 1i * rand(1,(N+1)^2 - 1);
            
            rwf = fields.RadiatingWaveField(k,alpha_n,beta_n);
            l2t = rwf.obtainFarField( 25, y);
                                              
            theta = pi * rand(100,1);
            phi = -pi + 2*pi*rand(100,1); 
            
            [F1,F2,F3] = l2t.eval( theta, phi);
            [V1, V2, V3] = rwf.evalFarFieldPattern( theta, phi, y );
                        
            err = sqrt( abs(F1 - V1).^2 + abs(F2 - V2).^2 + abs(F3 - V3).^2 );
            fprintf('Error: %e\n\n',max(err));            
            
            testcase.verifyLessThan(err,1e-5);
            
        end
        
    end
    
end