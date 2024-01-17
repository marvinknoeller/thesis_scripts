function [J,Curv_ints,Seg_ints] = Jacobian_FarField_maxwell_splineNonCirc(Var,p,Z,M,ts,R,S,T,aa,bb)
%Jacobian_FarField computes the Jacobian of the far field perturbation
%% Some help variables
num_x = length(p);
adder = zeros(3,3*num_x,num_x);
for kk=1:num_x
    adder(1:3,3*(kk-1)+1:3*kk,kk) = eye(3);
end
%% The Jacobian with respect to x1
% loop over the coordinates
parfor kk = 1:num_x
    for j = 1:3
        forParloop(:,:,kk,j) = derive_farfieldNonCirc(p,squeeze(adder(:,3*(kk-1)+j,:)),Var,M,Z,ts,R,S,T,aa,bb);
        [ints1(:,:,j,kk)] = Jacobian_Curv_Pen(p,squeeze(adder(:,3*(kk-1)+j,:)),M,num_x);
        [ints2(:,:,j,kk)] = Jacobian_Seg_Pen_New(p,squeeze(adder(:,3*(kk-1)+j,:)),M,num_x);
    end
end
%% 
% We have 180 points on the sphere. For different sizes this has to be
% changed.
J = zeros(180,3,3+3*(num_x-1));
Curv_ints = zeros(3,(num_x-1)*1*(M-1)+1,3*(num_x));
Seg_ints = zeros(num_x-1,3+3*(num_x-1));
for kk = 1:num_x
    for j =1:3
        J(:,:,j+3*(kk-1)) = forParloop(:,:,kk,j);
        Curv_ints(:,:,j+3*(kk-1)) = ints1(:,:,j,kk);
        Seg_ints(:,j+3*(kk-1)) = ints2(:,:,j,kk);
    end
end
