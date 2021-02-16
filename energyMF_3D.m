function [energyall,wfall]=energyMF_3D(ave,param)
%wfall(kz,NxNy,:)

sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
N=param.N;
t=param.t;
mu=param.mu;
energyall=cell(N(3),1);
wfall=cell(N(3),1);

Bx=spdiags([ones(N(1),1),ones(N(1),1)],[-1,1],N(1),N(1));   %banded mat for Nx
Bx(1,end)=1;Bx(end,1)=1;
By=spdiags([ones(N(2),1),ones(N(2),1)],[-1,1],N(2),N(2));   %banded mat for Ny
By(1,end)=1;By(end,1)=1;

Dmat=param.g*spdiags(ave,0,N(1)*N(2),N(1)*N(2));
zero=zeros(N(1)*N(2));

D_bdg=[zero,zero,zero,-Dmat;
         zero,zero,Dmat,zero;
         zero,conj(Dmat),zero,zero;
         -conj(Dmat),zero,zero,zero];
Z_bdg=kron(kron(sz,sz),speye(N(1)*N(2)))*param.Ez;


for kindex=1:N(3)
    k=estimate_k(param.klist(kindex),param);
    
    if k~=0
    K=-kron(Bx,speye(N(2)))*t(1)-kron(speye(N(1)),By)*t(2)+(2*t(1)+2*t(2)+param.energylist(kindex)-mu)*speye(N(1)*N(2));
    
    K_bdg=[K,zero,zero,zero;
           zero,K,zero,zero;
           zero,zero,-K,zero;
           zero,zero,zero,-K];
   
    H_bdg=K_bdg+D_bdg+Z_bdg;
    
   

    [vec,val]=eigs(H_bdg,k+10,'sm');
    val=real(diag(val));
    fprintf('kindex=%d,min(val)=%f,max(val)=%f\n',kindex,min(abs(val)),max(abs(val)));

    
    %restart with larger k   
    while max(val)<param.ED
        fprintf('k (%d) is too small, restart with k (%d)\n',k,k*2);
        k=k*2;
        [vec,val]=eigs(H_bdg,k,'sm');
        val=real(diag(val));
    end
    
    debyeindex=(abs(val)<=param.ED);
    val=val(debyeindex);
    vec=vec(:,debyeindex);
    [val,I]=sort(val);
    vec=vec(:,I);
    
    energyall{kindex}=val;
    wfall{kindex}=vec;
    end
end
end





