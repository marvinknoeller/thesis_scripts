function [D1, D2, F1] = computePenetrabDet(kR, eps_r)

kR_ins = sqrt(eps_r) * kR;

N=50;
D1 = zeros(1,N);
D2 = zeros(1,N);
F1 = zeros(1,N);

for n=1:N  
    
    j_n_out = sqrt(pi / (2*kR) ) * besselj(n+1/2,kR);
    j_n_min_1_out = sqrt(pi / (2*kR) ) .* besselj(n-1/2,kR);
    J_n_out = kR .* j_n_min_1_out - n * j_n_out;
    
    h_n = sqrt(pi / (2*kR) ) * besselh(n+1/2,1,kR);
    h_n_min_1 = sqrt(pi / (2*kR) ) .* besselh(n-1/2,1,kR);
    H_n = kR .* h_n_min_1 - n * h_n;
    
    j_n_in=  sqrt(pi / (2*kR_ins) ) * besselj(n+1/2,kR_ins);
    j_n_min_1_in = sqrt(pi / (2*kR_ins) ) .* besselj(n-1/2,kR_ins);
    J_n_in = kR_ins .* j_n_min_1_in - n * j_n_in;
    
    D1(n) = H_n .* j_n_in - h_n .* J_n_in;
    D2(n) = H_n .* j_n_in - 1/eps_r * h_n .* J_n_in;
    
    F1(n) = j_n_out .* J_n_in - J_n_out .* j_n_in; % ./ D1_n;
%     beta(index) = Einc.beta(index) .* ( 1/eps_r * j_n_out .* J_n_in - J_n_out .* j_n_in) ./ D2_n;
%     
%     gamma(index) = Einc.alpha(index) .* (j_n_out .* H_n - J_n_out .* h_n) ./ D1_n;
%     delta(index) = Einc.beta(index) .* ( j_n_out .* H_n - J_n_out .* h_n) ./ D2_n / sqrt(eps_r);
    
    
end

end