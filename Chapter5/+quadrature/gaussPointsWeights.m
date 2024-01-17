%> @file gaussPointsWeights.m
%> @brief The gaussPointsWeights() function.
% ======================================================================
%> @brief Computation of a Gaussian quadrature if the 3-term recurrence relation
%> for the orthogonal polynomials is known.
%>
%> This method is described in 
%>
%> Gene Golub and John H. Welsch, Mathematics of Computation 23 106 (1969), 221–230.
% ======================================================================
function [ x_n, w_n ] = gaussPointsWeights( alpha_n, beta_n, int_w )
%
% function [ x_n, w_n ] = gaussPointsWeights( alpha_n, beta_n, int_w )
%

J = diag(alpha_n) + diag(beta_n,1) + diag(beta_n,-1);

[V,D] = eig(J);

x_n = diag(D).';
w_n = abs(V(1,:)).^2 * int_w;

end

