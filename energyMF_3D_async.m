function [energyall,wfall]=energyMF_3D_async(ave,param)
%wfall(kz,NxNy,:)
tic;
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
pp.sx=sx;
pp.sy=sy;
pp.sz=sz;
pp.s0=s0;

N=param.N;
t=param.t;
mu=param.mu;
pp.param=param;

energyall=cell(N(3),1);
wfall=cell(N(3),1);

Bx=spdiags([ones(N(1),1),ones(N(1),1)],[-1,1],N(1),N(1));   %banded mat for Nx
Bx(1,end)=1;Bx(end,1)=1;
By=spdiags([ones(N(2),1),ones(N(2),1)],[-1,1],N(2),N(2));   %banded mat for Ny
By(1,end)=1;By(end,1)=1;

pp.Bx=Bx;
pp.By=By;

Dmat=param.g*spdiags(ave,0,N(1)*N(2),N(1)*N(2));
zero=zeros(N(1)*N(2));
pp.zero=zero;
D_bdg=[zero,zero,zero,-Dmat;
         zero,zero,Dmat,zero;
         zero,conj(Dmat),zero,zero;
         -conj(Dmat),zero,zero,zero];
Z_bdg=kron(kron(sz,sz),speye(N(1)*N(2)))*param.Ez;

pp.D_bdg=D_bdg;
pp.Z_bdg=Z_bdg;

kreq=arrayfun(@(x) estimate_k(x,param),param.klist(1:end/2));
pp.kreq=kreq;
knzlist=find(kreq);

% energyallknz=cell(length(knzlist),1);
% wfallknz=cell(length(knzlist),1);


f(1:length(knzlist)) = parallel.FevalFuture;
for idx = 1:length(knzlist)
    f(idx) = parfeval(@(x)diagonize(x,pp),3,knzlist(idx));
end

for idx = 1:length(knzlist)
    [completedIdx,val,vec,string] = fetchNext(f);
    energyall{knzlist(completedIdx)} = val;
    wfall{knzlist(completedIdx)} = vec;
    energyall{end-knzlist(completedIdx)}=val;
    wfall{end-knzlist(completedIdx)}=val;
    fprintf(string);
%     fprintf('Got result with index: %d.\n', completedIdx);
end

toc;
end

function [val,vec,string]=diagonize(kindex,pp)

    k=pp.kreq(kindex);
    K=-kron(pp.Bx,speye(pp.param.N(2)))*pp.param.t(1)-kron(speye(pp.param.N(1)),pp.By)*pp.param.t(2)...
        +(2*pp.param.t(1)+2*pp.param.t(2)+pp.param.energylist(kindex)-pp.param.mu)*speye(pp.param.N(1)*pp.param.N(2));
    
    K_bdg=[K,pp.zero,pp.zero,pp.zero;
           pp.zero,K,pp.zero,pp.zero;
           pp.zero,pp.zero,-K,pp.zero;
           pp.zero,pp.zero,pp.zero,-K];
   
    H_bdg=K_bdg+pp.D_bdg+pp.Z_bdg;
   

    [vec,val]=eigs(H_bdg,k+10,'sm');
    val=real(diag(val));
    s1=sprintf('kindex=%d,min(val)=%f,max(val)=%f,',kindex,min(abs(val)),max(abs(val)));

    
    %restart with larger k   
    while max(val)<pp.param.ED
        fprintf('k (%d) is too small, restart with k (%d)\n',k,k*2);
        k=k*2;
        [vec,val]=eigs(H_bdg,k,'sm');
        val=real(diag(val));
    end
    
    debyeindex=(abs(val)<=pp.param.ED);
    val=val(debyeindex);
    vec=vec(:,debyeindex);
    [val,I]=sort(val);
    vec=vec(:,I);
    s2=fprintf("len(val)=%d\n",length(val));
    string=strcat(s1,s2);
    
end




