%wfall(k,level,:)
function [energyall,wfall]=energyMF_r(ave,param)
Nk=param.Nk;
% energyall=zeros(Nk*4,1);
% wfall=zeros(Nk*4);
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;

B=spdiags([ones(Nk,1),ones(Nk,1)],[-1,1],Nk,Nk);
B(1,Nk)=1;
B(Nk,1)=1;

Dmat=spdiags(ave,0,Nk,Nk);
zero=zeros(Nk);
DD=[zero,zero,zero,Dmat;
 zero,zero,-Dmat,zero;
 zero,-conj(Dmat),zero,zero;
 conj(Dmat),zero,zero,zero];

T=kron(kron(sz,s0),B*(-param.t)+(2*param.t-param.mu)*speye(Nk));
Z=kron(kron(sz,sx),speye(Nk))*param.Ez;

D=DD*param.g/(2);

H=(T+D+Z);

[vec,val]=eig(full(H));
val=real(diag(val));
[val,I]=sort(val);
vec=vec(:,I);

energyall=val;
wfall=vec;
end


