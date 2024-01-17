function [der_R] = derive_par_wr_radius(t)
% Input: h is the current heigth
% Output: der_h is the derivative with respect to h, which is in my special
% case [0; 0; 1]

der_R = [cos(4*pi*t); sin(4*pi*t); zeros(1,length(t)) ];

% der_h = [zeros(1,length(t)); zeros(1,length(t)); ones(1,length(t)) ];
end
