%> @file gaussLegendre.m
%> @brief The gaussLegendre() function.
% ======================================================================
%> @brief Computation of the points and weights of a Gauss-Legendre
%> quadrature rule for a general interval [a,b]
%>
%> The function takes arguments @a N, @a a, @a b, where @a a, @a b are the
%> interval end points and a @a N+1 point rule is computed. The function
%> returns row vectors @a x and @a w of equal length containing the
%> quadrature points and weights, respectively.
%>
%> The function sets the coefficients of the 3 term recurrence relations of
%> the Legendre polynomials and uses the Golub/Welsh method as implemented
%> in gaussPointsWeights() to compute points and weights. Theses are then
%> scaled for the interval [a,b].
% ======================================================================
function [x, w] = gaussLegendre(N,a,b)
%
% function [x, w] = gaussLegendreQuad(N,a,b)
%

% Compute coefficients of 3 term recurrence relation
j = 1:N;
alpha_j = zeros(1,N+1);
beta_j = j ./ sqrt( (2*j-1) .* (2*j+1) );

% Compute points/weights for interval [-1,1]
[ x, w ] = quadrature.gaussPointsWeights( alpha_j, beta_j, 2 );

% Scale for [a,b]
x = a + (b-a)/2*(1+x);
w = (b-a)/2*w;
