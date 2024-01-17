function [chi_rel,HS,Cint,sigma_p,sigma_m] = chiral(F)

N = size(F,1)/2;
M = 1 / sqrt(2) * [eye(N), eye(N); 1i * eye(N), -1i * eye(N) ];
MM = inv(M);
FF = MM * F * M; 

Fpp = FF(1:N,1:N); 
Fmm = FF((N+1):2*N,(N+1):2*N);
Fpm = FF(1:N,(N+1):2*N);
Fmp = FF((N+1):2*N,1:N);
pp = svd(Fpp);
mm = svd(Fmm);
pm = svd(Fpm);
mp = svd(Fmp);
p = svd(FF);


chi = sqrt(norm( pp-mm,2)^2 + norm(pm - mp)^2);
Cint = norm(p,2)^2;
chi_rel = (chi) / sqrt(Cint);
HS = (norm(p,2)^2-2*norm(pp)*norm(mm)-2*norm(pm)*norm(mp))^(1/2)/sqrt(Cint);

sigma_p = sqrt(norm(pp)^2 + norm(mp)^2);
sigma_m = sqrt(norm(mm)^2 + norm(pm)^2);

end