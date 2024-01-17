%> @file sphericalHarmonics.m
%> @brief The sphericalHarmonics function.
% ======================================================================
%> @brief Evaluation of spherical harmonics
%>
%> The function evaluates spherical harmonics
%>
%> @f[
%>  Y_n^m(\hat{x}) = \sqrt{ \frac{ (2n+1) \, ( n - |m| )! }{ 4 \pi \, (n + |m|)! } } 
%>  \, P_n^{|m|}( \cos \vartheta ) \, \E^{\I m \varphi} 
%>  \, , \qquad n \in \NN_0 \, , \quad m = -n,\ldots,n \, .
%> @f]
%
% The function takes 3 input arguments @a theta, @a phi, @a n. Here, @a theta
% and @a phi are row vectors of equal length specifying the evaluation
% points on the unit sphere by
%>
%> @f[
%>  \hat{x} = ( \cos(\varphi) \, \sin(\vartheta) \, ,  \, (\sin(\varphi) \,
%>  \sin(\vartheta) \, , \, \cos(\vartheta) )^\top
%> @f]
%>
%> The argument @n specifies the degree of the spherical harmonics to
% evaluate. The function returns a matrix of size 2n+1 x length(theta) matrix 
% containing the of values
%  [ Y_n^{-n}(theta,phi); Y_n^{-n+1}(theta,phi); ... Y_n^n(theta,phi) ]
% ======================================================================
function Y = sphericalHarmonics( theta, phi, n)
%
% function Y = sphericalHarmonics( theta, phi, n)
%

m = (-n:n).';

assocLegendre = 1/sqrt(2*pi) * legendre(n, cos(theta), 'norm');
assocLegendre = [ flip(assocLegendre,1); assocLegendre(2:end,:) ];
% assocLegendre = [(-1).^(1:n).* flip(assocLegendre,1); assocLegendre(2:end,:) ];
Y =  exp(1i * m * phi) .* assocLegendre;
% Y =  (-1).^(m(2:end)) .* exp(1i * m * phi) .* assocLegendre;

end