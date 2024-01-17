%> @file VectorField.m
%> @brief Contains the VectorField class.
% ======================================================================
%> @brief Abstract class representing a vector field that can be evaluated.
% ======================================================================
classdef VectorField < handle
    % classdef VectorField < handle
    
    methods (Abstract)
        
        %> @brief Evaluation of the vector field
        %>
        %> Returns matrices V1, V2, V3 of the same size as X1, X2, X3
        %> containing the values of the components of the vector fields in
        %> the points given by the Xj
        %>
        %> @param[in] vf The VectorField instance.
        %> @param[in] X1 Matrix of x1 coordinates of the points at which
        %>   evaluation takes place.
        %> @param[in] X2 Matrix of x2 coordinates, same size as X1.
        %> @param[in] X3 Matrix of x3 coordinates, same size as X1.
        %> @param[out] V1 First component of result, same size as X1.
        %> @param[out] V2 Second component of result, same size as X1.
        %> @param[out] V3 Third component of result, same size as X1.
        [V1, V2, V3] = eval( vf, X1, X2, X3 );
        
    end
    
end

