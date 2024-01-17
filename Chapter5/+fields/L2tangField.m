%> @file L2tangField.m
%> @brief Contains the L2tangField class.
% ======================================================================
%> @brief Class representing a square integrable tangential vector field on
%> the unit sphere.
%>
%> A tangential field on the unit sphere can be defined as a superposition
%> of vector spherical harmonics,
%>
%> @f[
%> W(\hat{x}) = \sum_{n=1}^N \sum_{m = -n}^n \left[ a_n^m \, U_n^m (\hat{x}) 
%> + b_n^m \, V_n^m (\hat{x}) ]
%> @f]
%>
%> The tangential field is thus defined
%> by its coefficients a_n^m and b_n^m
% ======================================================================
classdef L2tangField < handle
    % classdef L^2tangField < handle
    
    properties (SetAccess = protected)
        
        %> @brief The degree of the highest vector spherical harmonic in
        %> the expansion.
        degree = 1;
        
        %> @brief The coefficients for the functions @f$ U_n^m @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like
        %>
        %> a_1^{-1} a_1^{0} a_1^{1} a_2^{-2} ... a_2^{2} a_3^{-3} ... a_degree^{degree}
        %>
        a = [0, 0, 0];
        
        %> @brief The coefficients for the functions @f$ V_n^m @f$.
        %>
        %> These form a row vector of length degree^2 - 1 and are ordered
        %> like the coefficients in @a a.
        b = [0, 0, 0];
        
    end
    
    methods
        
        %> @brief Construction
        %>
        %> A call with less than two arguments initializes the object to
        % the zero function. If two argumnts @a a, @b a are provided, these
        % are define the basis coefficients. They should be row vectors of 
        % equal length n^2-1. The number n  then is the degree.
        function l2t = L2tangField( a, b )
            
            if ( nargin >= 2 )
                l2t.setCoeffs( a, b );
            end
            
        end
        
        %> @brief Set degree and expanson coefficients
        function setCoeffs(l2t, a, b)
            
            %> @todo Check sizes
            
            l2t.a = a;
            l2t.b = b;
            l2t.degree = sqrt( length(a) + 1 ) - 1;            
            
        end
        
        %> @brief Evaluation of the function.
        %>
        %> The function computes the evaluation of the L2tangField in the
        %> points
        %>
        %> @f[
        %>  \hat{x} = ( \cos(\varphi) \, \sin(\vartheta) \, ,  \, (\sin(\varphi) \,
        %>  \sin(\vartheta) \, , \, \cos(\vartheta) )^\top
        %> @f]
        %> 
        %> on the unit sphere.
        %>
        %> The arguments Theta and Phi should be arrays of equal size with
        %> values in [0,\pi] and [-pi,pi], respectively.
        %>
        %> The function returns arrays @a V1, @a  V2, @a V3 of the same
        %> size as @a Theta and @a Phi which contain the components of the
        %> field in the evaluation points.
        function [V1, V2, V3] = eval( l2t, Theta, Phi )
            
            theta = Theta(:).';
            phi = Phi(:).';

            V1 = zeros(1,length(theta));
            V2 = zeros(1,length(theta));
            V3 = zeros(1,length(theta));

            for n=1:l2t.degree
                [U, V] = specialfunc.vectorSphericalHarmonics( theta, phi, n );
                
                mRange = n^2:((n+1)^2-1);
                
                V1 = V1 + ( l2t.a(mRange) * U(:,:,1) ) + ( l2t.b(mRange) * V(:,:,1) );
                V2 = V2 + ( l2t.a(mRange) * U(:,:,2) ) + ( l2t.b(mRange) * V(:,:,2) );
                V3 = V3 + ( l2t.a(mRange) * U(:,:,3) ) + ( l2t.b(mRange) * V(:,:,3) );

            end
            
            V1 = reshape( V1, size(Theta) );
            V2 = reshape( V2, size(Theta) );
            V3 = reshape( V3, size(Theta) );

        end
    
    end
    
    methods (Static)
        
        %> @brief Compute the best approximation on the unit sphere of a given 
        %> VectorField by a L2tangField.
        %>
        %> The approximation is done by orthogonal projection onto
        %> the space spanned by the vector spherical harmonics of degree less than or
        %> equal to @a N in the L^2-norm.
        function ewf = bestApprox(vf, N)
        
            % obtain the quadrature points and weights
            N_theta = 50;
            N_phi = 60;
            
            [Theta, Phi, W] = quadrature.sphereGaussLegTrap( N_theta, N_phi );            
            
            % evaluate the VectorField in the quadrature points
            if ( isa(vf,'fields.VectorFieldOnSphere') )
                [ F1, F2, F3 ] = vf.eval(Theta,Phi);
            else                
                [ F1, F2, F3 ] = vf.eval( ...
                    cos(Phi) .* sin(Theta), ...
                    sin(Phi) .* sin(Theta), ...
                    cos(Theta) );
            end
            
            % intialize coefficients
            alpha = zeros(1,(N+1)^2 - 1);
            beta = zeros( size(alpha) );
            
            % for each vector spherical harmonic up to degree N
            for n=1:N
                
                Wmat = W * ones(1,2*n+1);
                
                % evaluate the vsh on the quadrature points
                [U, V] = specialfunc.vectorSphericalHarmonics( Theta, Phi, n);
                
                % scalar producs in L^2_tangential
                index = n^2:( (n+1)^2 - 1 );
                beta(index)  = F1 * ( V(:,:,1)' .* Wmat ) + F2 * ( V(:,:,2)' .* Wmat ) + F3 * ( V(:,:,3)' .* Wmat );
                alpha(index) = F1 * ( U(:,:,1)' .* Wmat ) + F2 * ( U(:,:,2)' .* Wmat ) + F3 * ( U(:,:,3)' .* Wmat );
                
            end
            
            % define the EntireWaveField
            ewf = fields.L2tangField;
            ewf.setCoeffs( alpha, beta );
            
        end
        
    end
    
end