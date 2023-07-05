function [S] = AdaptiveStructureGraph (X_lib,X,kmax)
kmax = kmax+1;
kmin = round(kmax/10)+1;
[idx, distX] = knnsearch(X_lib,X,'k',kmax);
[N,~] = size(X);
degree_x = tabulate(idx(:));
kmat = degree_x(:,2);
kmat(kmat >= kmax)=kmax;
kmat(kmat <= kmin)=kmin;
if length(kmat) < N
    kmat(length(kmat)+1:N) = kmin;
end
S = zeros(N,N);
for i = 1:N
    K = kmat(i);
    k = K-1;
    id_x = idx(i,1:K);
    di = distX(i,1:K);
    W = (di(K)-di)/(k*di(K)-sum(di(1:k))+eps);
    S(i,id_x) = W;
end
