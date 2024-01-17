%> @file WrapperFieldOnSphere.m
%> @brief Contains the WrapperFieldOnSphere class.
% ======================================================================
%> @brief Wrapper class that passes evaluation to another function.
% ======================================================================
classdef WrapperFieldOnSphere < fields.VectorFieldOnSphere
   
    
    properties (SetAccess = protected)
        
        evalFunc = [];
        
    end
    
    methods
        
        %> @brief Construction
        %>
        %> If at least one argument, the first argument is a handle to the
        %> function the eval() call is supposed to be passed to.
        function wf = WrapperFieldOnSphere( evalHandle )
            
            if ( nargin >= 1 )
                wf.evalFunc = evalHandle;
            end
            
        end
        
        %> @brief Set center of expansion 
        function setEvalFunction (wf, evalHandle)
            
            wf.evalFunc = evalHandle;
            
        end
        
        
        %> @brief Evaluate by passing call to handle function
        function [V1, V2, V3] = eval( wf, Theta, Phi )
            
            [V1, V2, V3] = wf.evalFunc( Theta, Phi );
            
        end
        
    end
    
end