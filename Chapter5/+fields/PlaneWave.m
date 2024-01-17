%> @file PlaneWave.m
%> @brief Implementation of a VectorField as a plane wave.
% ======================================================================
%> @brief Simple implementation of VectorField as a plane wave.
%>
%> A plane wave is defined by wave number, direction of propagation and
%> polarization. The last two are column vectors of length three and should
%> be perpendicular.
% ======================================================================
classdef PlaneWave < fields.VectorField
    
    
    properties (SetAccess = protected)
        
        %> @brief The wave number
        wavenum = 1;
        
        %> @brief Direction of propagation
        d = [1;0;0];
        
        %> @brief Polarization
        p = [0;1;0];
               
    end
    
    methods
        
        %> @brief Construction with optional initialization of parameters.
        %
        % The constructor takes arguments @a wavenum (the wave number), 
        % @a d (direction of incidence), @a p (polarization vector). The
        % the two latter arguments should be column vectors of length 3.
        % The vector d must have Euclidean norm 1 and @a d and @a p must be
        % orthogonal. Default values, if these parameters are omited, are 
        % @a wavenum = 1, @a d = [1;0;0] and @a p = [0;1;0].
        function pw = PlaneWave( wavenum, d, p )
            
            if ( nargin > 0 )
                pw.wavenum = wavenum;
            end
            
            if ( nargin > 1 )
                pw.d = d;
            end
            
            if ( nargin > 2 )
                pw.p = p;
            end
           
            if ( abs( norm(pw.d) - 1.0 ) > 1e-14 )
                error('PlaneWave construction: direction of propagation must have norm 1.');
            end
           
            if ( abs( dot(pw.d,pw.p) ) > 1e-14 )
                error('PlaneWave construction: direction of propagation and polarization must be orthogonal.');
            end
            
        end
        
        %> @brief Implementation of VectorField.eval()
        function [V1, V2, V3] = eval( pw, X1, X2, X3 )
            
            exponential = exp( 1i * pw.wavenum * ( pw.d(1) * X1 + pw.d(2) * X2 + pw.d(3) * X3 ) );
            
            V1 = pw.p(1) * exponential;
            V2 = pw.p(2) * exponential;
            V3 = pw.p(3) * exponential;
            
        end
        
    end
    
end