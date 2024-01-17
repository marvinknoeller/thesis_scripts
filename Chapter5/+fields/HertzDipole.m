%> @file HertzDipole.m
%> @brief Contains the HertzDipole class.
% ======================================================================
%> @brief Implementation of a VectorField as a Hertz dipole.
%>
%> The HertzDipole represents the electric field due to a Hertz dipole
%> with position y and polarization p.
% ======================================================================
classdef HertzDipole < fields.VectorField
    
    
    properties (SetAccess = protected)
        
        %> @brief The wave number
        wavenum = 1;
        
        %> @brief Source point
        y = [0;0;0];
        
        %> @brief Polarization vector
        p = [0;0;1];
               
    end
    
    methods
        
        function hd = HertzDipole( wavenum, y, p )
            
            if ( nargin > 0 )
                hd.wavenum = wavenum;
            end
            
            if ( nargin > 1 )
                hd.y = y;
            end
            
            if ( nargin > 2 )
                hd.p = p;
            end
            
        end
        
        function [V1, V2, V3] = eval( hd, X1, X2, X3 )
            
            Z1 = X1 - hd.y(1);
            Z2 = X2 - hd.y(2);
            Z3 = X3 - hd.y(3);
            
            % cross product Z x p
            c1 = Z2 * hd.p(3) - Z3 * hd.p(2);
            c2 = Z3 * hd.p(1) - Z1 * hd.p(3);
            c3 = Z1 * hd.p(2) - Z2 * hd.p(1);
            
            % cross product (Z x p) x Z
            d1 = c2 .* Z3 - c3 .* Z2;
            d2 = c3 .* Z1 - c1 .* Z3;
            d3 = c1 .* Z2 - c2 .* Z1;
            
            % vector p - 3 * (P*hZ) hZ
            R = sqrt( Z1.^2 + Z2.^2 + Z3.^2 );
            sp = Z1*hd.p(1) + Z2*hd.p(2) + Z3*hd.p(3);
            e1 = hd.p(1) - 3 * sp .* Z1 ./ R.^2;
            e2 = hd.p(2) - 3 * sp .* Z2 ./ R.^2;
            e3 = hd.p(3) - 3 * sp .* Z3 ./ R.^2;
            
            eFactor = (1i * hd.wavenum * R - 1) / hd.wavenum^2;
            exponential = 1i * hd.wavenum * exp( 1i * hd.wavenum * R ) ./ (4*pi * R.^3);
            
            V1 = exponential .* ( d1 + eFactor .* e1 );
            V2 = exponential .* ( d2 + eFactor .* e2 );
            V3 = exponential .* ( d3 + eFactor .* e3 );
            
        end
        
        
        function [hdinf1, hdinf2, hdinf3] = farFieldPattern( hd, Theta, Phi )
            
            theta = Theta(:);
            phi = Phi(:);
            
            Z1 = cos(phi) .* sin(theta);
            Z2 = sin(phi) .* sin(theta);
            Z3 = cos(theta);
            
            % cross product Z x p
            c1 = Z2 * hd.p(3) - Z3 * hd.p(2);
            c2 = Z3 * hd.p(1) - Z1 * hd.p(3);
            c3 = Z1 * hd.p(2) - Z2 * hd.p(1);
            
            % cross product (Z x p) x Z
            d1 = c2 .* Z3 - c3 .* Z2;
            d2 = c3 .* Z1 - c1 .* Z3;
            d3 = c1 .* Z2 - c2 .* Z1;
            
            phaseFactor = exp( -1i * hd.wavenum * ( hd.y(1)*Z1 + hd.y(2)*Z2 + hd.y(3)*Z3 ) );
            
            hdinf1 = 1i * hd.wavenum * phaseFactor .* d1;
            hdinf2 = 1i * hd.wavenum * phaseFactor .* d2;
            hdinf3 = 1i * hd.wavenum * phaseFactor .* d3;
            
            hdinf1 = reshape( hdinf1, size(Theta) );
            hdinf2 = reshape( hdinf2, size(Theta) );
            hdinf3 = reshape( hdinf3, size(Theta) );
            
        end
        
    end
    
end