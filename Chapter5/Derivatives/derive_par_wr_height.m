function [der_h] = derive_par_wr_height(t)
% Input: h is the current heigth
% Output: der_h is the derivative with respect to h, which is in my special
% case [0; 0; 1]

der_h = [zeros(1,length(t)); zeros(1,length(t)); t./1 - 1/2 ];

% der_h = [zeros(1,length(t)); zeros(1,length(t)); ones(1,length(t)) ];
end

