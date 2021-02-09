%wfall(k,level,:)
function [energyall,wfall]=energyMF(ave,param)
Nk=param.Nk;
energyall=zeros(Nk,4);
wfall=zeros(Nk,4,4);
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
T=zeros(4,4,Nk);
D=zeros(4,4,Nk);
Z=zeros(4,4,Nk);

for kindex=1:Nk
    T(:,:,kindex)=kron(sz,s0)*param.energylist(kindex);
    D(:,:,kindex)=kron(sy,sy)*ave*param.g/(2*param.Nk);
    Z(:,:,kindex)=param.Ez*kron(sz,sx);
end

H=T+D+Z;

for kindex=1:Nk
    [vec,val]=eig(H(:,:,kindex));
    val=real(diag(val));
    [val,I]=sort(val);
    vec=vec(:,I);
    energyall(kindex,:)=val;
    for ii=1:4
        wfall(kindex,ii,:)=vec(:,ii);
    end
end

