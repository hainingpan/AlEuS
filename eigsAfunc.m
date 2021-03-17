% [L,U,p]=lu(H_bdg,'vector');
dA=decomposition(H_bdg,'lu');
Afun=@(x) dA\x;
% Afun=@(x) U\(L\(x(p)));
[vvec,vval]=eigs(Afun,length(H_bdg),k/2,'lr','Display',1,'IsFunctionSymmetric',1,'SubspaceDimension',2*k);
vval=1./diag(vval);