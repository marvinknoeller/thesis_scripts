%> @file vectorSphericalHarmonics.m
%> @brief The vectorSphericalHarmonics function.
% ======================================================================
%> @brief Evaluation of vector spherical harmonics
%>
%> The function evaluates vector spherical harmonics
%>
%> @f[
%>  U_n^m(\hat{x}) = \frac{1}{\sqrt{ n \, (n+1) } } \, \operatorname{Grad}_{S^2}
%>  Y_n^m(\hat{x}) \, ,
%>  \qquad 
%>  V_n^m(\hat{x}) = \hat{x} \times U_n^m(\hat{x})
%> @f]
%>
%> in pointts @f$ \hat{x} @f$ on the unit sphere given by parameters phi,
%> theta with
%>
%> @f[
%>  \hat{x} = ( \cos(\varphi) \, \sin(\vartheta) \, ,  \, (\sin(\varphi) \,
%>  \sin(\vartheta) \, , \, \cos(\vartheta) )^\top
%> @f]
%>
%> The function takes parameters @a theta, @a phi which are row vectors
%> of equal length and a natural number @a n as input. The ouput are
%> two arrays @a U, @a V of dimension 3. The first dimension has length 2n+1
%> and corresponds to the parameter m, the second dimension has length equal
%> to length(theta) and corresponds to the evaluation point, the third
%> dimension has length 3 and corresponds to the component of the vector
%> field.
% ======================================================================
function [U, V] = vectorSphericalHarmonics( theta, phi, n )

assocLegendre = zeros(n+2,length(theta));

assocLegendre(1:n+1,:) = legendre(n, cos(theta), 'norm');

m = (0:n).';
exp_term = exp(1i * m * phi) / sqrt( 2*pi * n * (n+1) );

coeff_theta = zeros(2*n+1,length(theta));
coeff_phi = zeros(2*n+1,length(theta));

coeff_theta(n+1:end,:) = exp_term .* ( (m * cot(theta)) .* assocLegendre(1:n+1,:) ...
    - ( sqrt( ( n-m ) .* ( n+m+1 ) ) * ones(size(theta)) ) .* assocLegendre(2:n+2,:) );
coeff_phi(n+1:end,:) = exp_term .* ( 1i * m * (1./ sin(theta)) ) .* assocLegendre(1:n+1,:);

index = ( sin(theta) < 1e-10 );
coeff_theta(:,index) = 0;
coeff_theta(n+2,index) = n*(n+1)/2 * sqrt( (2*n+1) / 2 / n / (n+1) ) * exp_term(2,index);
coeff_phi(:,index) = 0;
coeff_phi(n+2,index) = 1i * n*(n+1)/2 * sqrt( (2*n+1) / 2 / n / (n+1) ) * exp_term(2,index);
% coeff_phi(n+2,index) = 1i * m * n*(n+1)/2 * sqrt( (2*n+1) / 2 / n / (n+1) ) * exp_term(2,index);
index = abs(theta - pi) < 1e-10;
coeff_theta(n+2,index) = (-1).^n * coeff_theta(n+2,index);
coeff_phi(n+2,index) = (-1).^(n+1) * coeff_phi(n+2,index);


coeff_theta(1:n,:) = conj( flip(coeff_theta(n+2:end,:),1) );
coeff_phi(1:n,:) =  conj( flip(coeff_phi(n+2:end,:),1) );

% coeff_theta(1:n,:) = (-1).^(m(2:end)) .* conj( flip(coeff_theta(n+2:end,:),1) );
% coeff_phi(1:n,:) =  (-1).^(m(2:end)) .* conj( flip(coeff_phi(n+2:end,:),1) );

U = zeros(2*n+1,length(theta),3);
V = zeros(2*n+1,length(theta),3);

U(:,:,1) = coeff_theta .* ( ones(2*n+1,1) * ( cos(phi) .* cos(theta) ) ) ...
    - coeff_phi .* ( ones(2*n+1,1) * sin(phi) );
U(:,:,2) = coeff_theta .* ( ones(2*n+1,1) * ( sin(phi) .* cos(theta) ) ) ...
    + coeff_phi .* ( ones(2*n+1,1) * cos(phi) );
U(:,:,3) = coeff_theta .* ( ones(2*n+1,1) * ( -sin(theta) ) );

V(:,:,1) = ones(2*n+1,1) * ( sin(phi) .* sin(theta) ) .* U(:,:,3) ...
    - ones(2*n+1,1) * ( cos(theta) ) .* U(:,:,2);
V(:,:,2) = ones(2*n+1,1) * ( cos(theta) ) .* U(:,:,1) ...
    - ones(2*n+1,1) * ( cos(phi) .* sin(theta) ) .* U(:,:,3);
V(:,:,3) = ones(2*n+1,1) * ( cos(phi) .* sin(theta) ) .* U(:,:,2) ...
    - ones(2*n+1,1) * ( sin(phi) .* sin(theta) ) .* U(:,:,1);

end