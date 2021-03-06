function [energyall,wfall]=energyMF_3D_sync(ave,wfall,param)
%wfall(kz,NxNy,:)
tic;
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
N=param.N;
t=param.t;
mu=param.mu;
energyall=cell(N(3),1);
% wfall=cell(N(3),1);

Bx=spdiags([ones(N(1),1),ones(N(1),1)],[-1,1],N(1),N(1));   %banded mat for Nx
By=spdiags([ones(N(2),1),ones(N(2),1)],[-1,1],N(2),N(2));   %banded mat for Ny
if param.pbc==1
    Bx(1,end)=1;Bx(end,1)=1;
    By(1,end)=1;By(end,1)=1;
end

Dmat=param.g*spdiags(ave,0,N(1)*N(2),N(1)*N(2));
zero=zeros(N(1)*N(2));

D_bdg=[zero,zero,zero,-Dmat;
         zero,zero,Dmat,zero;
         zero,conj(Dmat),zero,zero;
         -conj(Dmat),zero,zero,zero];
Z_bdg=kron(kron(sz,sz),speye(N(1)*N(2)))*param.Ez;

kreq=arrayfun(@(x) estimate_k(x,param),param.klist(1:end/2));

knzlist=find(kreq);
% knzlist=1:length(kreq);

energyallknz=cell(length(knzlist),1);
wfallknz=cell(length(knzlist),1);
parfor knzindex=1:length(knzlist)
    kindex=knzlist(knzindex);
    k=kreq(knzlist(knzindex));
    K=-kron(Bx,speye(N(2)))*t(1)-kron(speye(N(1)),By)*t(2)+(2*t(1)+2*t(2)+param.energylist(kindex)-mu)*speye(N(1)*N(2))...
        +spdiags(param.muVar,0,N(1)*N(2),N(1)*N(2));
    
    K_bdg=[K,zero,zero,zero;
           zero,K,zero,zero;
           zero,zero,-K,zero;
           zero,zero,zero,-K];
   
    H_bdg=K_bdg+D_bdg+Z_bdg;
    
   
    H_bdg=real((H_bdg+H_bdg')/2);
    k=max(20,2*ceil(1.05*k/2));
    if isempty(wfall{kindex})
        [vec,val]=eigs(H_bdg,k,'sm');
    else
        [vec,val]=eigs(H_bdg,k,'sm','StartVector',wfall{kindex}(:,1));
    end
        
    val=real(diag(val));
    fprintf('kindex=%d,min(val)=%f,max(val)=%f,',kindex,min(abs(val)),max(abs(val)));

    
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
    fprintf("len(val)=%d\n",length(val));
    
    energyallknz{knzindex}=val;
    wfallknz{knzindex}=vec;

end

for knzindex=1:length(knzlist)
    
    energyall{knzlist(knzindex)}=energyallknz{knzindex};
    wfall{knzlist(knzindex)}=wfallknz{knzindex};
    energyall{end-knzlist(knzindex)+1}=energyall{knzlist(knzindex)};
    wfall{end-knzlist(knzindex)+1}=wfall{knzlist(knzindex)};
end
toc;
end





