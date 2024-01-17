function plotOnSphere(vf, s)
%
% Plot the VectorFielf vf on the Sphere s
%
%  6 Subplots are created, real and imaginary part of each cartesean 
% component of vf, respectively.
%

% obtain the quadrature points and weights
if (~ s.isQuadDefined )
    s.predefQuadPoints(100,100);
end

[X1, X2, X3] = getCartesianQuadratureCoord(s);

if ( isa(vf,'L2tangField') )
    [Theta, Phi] = s.getAngularQuadratureCoord;
    [V1, V2,V3] = vf.eval( Theta, Phi );
else
    [V1, V2, V3] = vf.evalOnSurf( s );
end

X1 = reshape(X1,2*s.Nazim, s.Npolar+1);
X2 = reshape(X2,2*s.Nazim, s.Npolar+1);
X3 = reshape(X3,2*s.Nazim, s.Npolar+1);
V1 = reshape(V1,2*s.Nazim, s.Npolar+1);
V2 = reshape(V2,2*s.Nazim, s.Npolar+1);
V3 = reshape(V3,2*s.Nazim, s.Npolar+1);

figure

subplot(2,3,1)
surf(X1,X2,X3,real(V1))
shading interp
axis equal
colorbar

subplot(2,3,4)
surf(X1,X2,X3,imag(V1))
shading interp
axis equal
colorbar

subplot(2,3,2)
surf(X1,X2,X3,real(V2))
shading interp
axis equal
colorbar

subplot(2,3,5)
surf(X1,X2,X3,imag(V2))
shading interp
axis equal
colorbar

subplot(2,3,3)
surf(X1,X2,X3,real(V3))
shading interp
axis equal
colorbar

subplot(2,3,6)
surf(X1,X2,X3,imag(V3))
shading interp
axis equal
colorbar
