
%> @file sphereGaussLegTrap.m
%> @brief The sphereGaussLegTrap() function.
% ======================================================================
%> @brief Compute quadrature points (in polar and azimuthal angles) on
% the unit sphere for quadrature by a tensor product quadrature rule, with 
%> a Gauss-Legendre rule in the polar angle and a composite trapezoidal rule 
%> in the azimuthal angle. 2*Nazim points in the azimuthal and Npolar points 
%> in the polar direction are used.
% ======================================================================
function [Theta, Phi, W] = sphereGaussLegTrap(Npolar, Nazim)

[ t, w_t ] = quadrature.gaussLegendre(Npolar-1,0,pi);
p = linspace(-pi,pi,2*Nazim+1);
p = p(2:end);

[Theta, Phi] = meshgrid(t,p);
Theta = Theta(:).';
Phi = Phi(:).';

W = pi / Nazim * ones(length(p),1) * ( w_t .* sin(t) );
W = W(:);

end

