%> @file VectorFieldOnSphere.m
%> @brief Contains the VectorFieldOnSphere class.
% ======================================================================
%> @brief Abstract class representing a vector field that can be evaluated
%> on the unit sphere.
% ======================================================================
classdef VectorFieldOnSphere < handle
    
    methods (Abstract)
        
        %> @brief Evaluation of the vector field in a point specified by
        %> polar coordinates.
        %>
        %> Returns matrices V1, V2, V3 of the same size as Theta, Phi
        %> containing the values of the components of the vector fields in
        %> the points given by the angulat coordinates Theta, Phi
        %>
        %> @param[in] vf The VectorFieldOnSphere instance.
        %> @param[in] Theta Matrix of polar angular coordinates of the points at which
        %>   evaluation takes place.
        %> @param[in] Phi Matrix of azimuthal angular coordinates, same size as Theta.
        %> @param[out] V1 First component of result, same size as Theta.
        %> @param[out] V2 Second component of result, same size as Theta.
        %> @param[out] V3 Third component of result, same size as Theta.
        [V1, V2, V3] = eval( vf, Theta, Phi );
        
    end
    
end

