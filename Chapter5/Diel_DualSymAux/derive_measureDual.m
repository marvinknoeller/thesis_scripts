function [derivative] = derive_measure(F,H)

N = size(F,1)/2;
M = 1 / sqrt(2) * [eye(N), eye(N); 1i * eye(N), -1i * eye(N) ];
MM = inv(M);
FF = MM * F * M;

Fpp = FF(1:N,1:N);
Fmm = FF((N+1):2*N,(N+1):2*N);
Fpm = FF(1:N,(N+1):2*N);
Fmp = FF((N+1):2*N,1:N);

HH = MM * H * M;


Hpp = HH(1:N,1:N);
Hmm = HH((N+1):2*N,(N+1):2*N);
Hpm = HH(1:N,(N+1):2*N);
Hmp = HH((N+1):2*N,1:N);

pp = svd(Fpp);
mm = svd(Fmm);
pm = svd(Fpm);
mp = svd(Fmp);
p = svd(FF);

h_F = norm(pp)*norm(mm) + norm(pm)*norm(mp);
zae2 = real(trace(Hpp'*Fpp))/norm(pp) * norm(mm) + real(trace(Hpm'*Fpm))/norm(pm) * norm(mp)...
    + real(trace(Hmp'*Fmp))/norm(mp) * norm(pm) + real(trace(Hmm'*Fmm))/norm(mm) * norm(pp);
chiHS = ((norm(p)^2 - 2*h_F))^(1/2);
derivative = ((real(trace(HH'*FF)) - zae2) * norm(p)^2 - real(trace(HH'*FF))*chiHS^2)/(chiHS*norm(p)^3);



end