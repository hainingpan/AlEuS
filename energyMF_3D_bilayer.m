function [energyall,wfall]=energyMF_3D_bilayer(ave,wfall,param)
%wfall(kz,NxNy,:)
% tic;
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
N_Al=param.N_Al;
N_FM=param.N_FM;
N=param.N;
t_Al=param.t_Al;
t_FM=param.t_FM;
t_int=param.t_int;
mu_Al=param.mu_Al;
mu_FM=param.mu_FM;

energyall=cell(N(3),1);
% wfall=cell(N(3),1);

%Hopping for Al
Bx_Al=spdiags([ones(N_Al(1),1),ones(N_Al(1),1)],[-1,1],N_Al(1),N_Al(1));   %banded mat for Nx_Al
By_Al=spdiags([ones(N_Al(2),1),ones(N_Al(2),1)],[-1,1],N_Al(2),N_Al(2));   %banded mat for Ny_Al
if param.pbc==1
    Bx_Al(1,end)=1;Bx_Al(end,1)=1;
    By_Al(1,end)=1;By_Al(end,1)=1;
end

%Hopping for FM
Bx_FM=spdiags([ones(N_FM(1),1),ones(N_FM(1),1)],[-1,1],N_FM(1),N_FM(1));   %banded mat for Nx_FM
By_FM=spdiags([ones(N_FM(2),1),ones(N_FM(2),1)],[-1,1],N_FM(2),N_FM(2));   %banded mat for Ny_FM
if param.pbc==1
    Bx_FM(1,end)=1;Bx_FM(end,1)=1;
    By_FM(1,end)=1;By_FM(end,1)=1;
end


%Hopping between Al and FM
if N_FM(1)>0
    Bx_int=sparse([N_FM(1)],[1],[1],N_FM(1),N_Al(1));
    By_int=speye(N_FM(2)); %Assume N_FM(2)==N_Al(2)
    idx_FM=speye(N_FM(1)); 
    idx_FM(end,end)=0;
    intx_FM=sparse(N_FM(1),N_FM(1)); 
    intx_FM(end,end)=1;
else
    Bx_int=sparse(N_FM(1),N_Al(1));
    By_int=speye(N_FM(2));
    idx_FM=speye(N_FM(1)); 
    intx_FM=sparse(N_FM(1),N_FM(1)); 
end

idx_Al=speye(N_Al(1));
idx_Al(1,1)=0;

intx_Al=sparse(N_Al(1),N_Al(1));
intx_Al(1,1)=1;

K_int=-kron(Bx_int,By_int)*t_int(1);

Dmat=param.g*spdiags(ave,0,N(1)*N(2),N(1)*N(2));
% zero=zeros(N(1)*N(2));


% D_bdg=[zero,zero,zero,-Dmat;
%          zero,zero,Dmat,zero;
%          zero,conj(Dmat),zero,zero;
%          -conj(Dmat),zero,zero,zero];
     
D_bdg=kron(sparse([1],[4],[1],4,4),-Dmat)...
    +kron(sparse([2],[3],[1],4,4),Dmat)...
    +kron(sparse([3],[2],[1],4,4),conj(Dmat))...
    +kron(sparse([4],[1],[1],4,4),-conj(Dmat));
     
Z_bdg=kron(kron(sz,sz),spdiags([ones(N_FM(1)*N_FM(2),1);zeros(N_Al(1)*N_Al(2),1)],0,N(1)*N(2),N(1)*N(2)))*param.Ez;

kreq_est=arrayfun(@(x) estimate_k_bilayer(x,param),param.klist(1:end/2));
kreq=max(kreq_est,20);


knzlist=find(kreq);
% knzlist=1:length(kreq);

energyallknz=cell(length(knzlist),1);
wfallknz=cell(length(knzlist),1);
for knzindex=1:length(knzlist)
    kindex=knzlist(knzindex);
    k=kreq(knzlist(knzindex));
    
    K_Al=-kron(Bx_Al,speye(N_Al(2)))*t_Al(1)-kron(speye(N_Al(1)),By_Al)*t_Al(2)...
        +(2*t_Al(1)+2*t_Al(2)+param.energylist_Al(kindex)-mu_Al)*kron(idx_Al,speye(N_Al(2)))...
        +(t_Al(1)+t_int(1)+2*t_Al(2)+param.energylist_Al(kindex)-mu_Al)*kron(intx_Al,speye(N_Al(2)))...
        +spdiags(param.muVar,0,N_Al(1)*N_Al(2),N_Al(1)*N_Al(2));
    
    K_FM=-kron(Bx_FM,speye(N_FM(2)))*t_FM(1)-kron(speye(N_FM(1)),By_FM)*t_FM(2)...
        +(2*t_FM(1)+2*t_FM(2)+param.energylist_FM(kindex)-mu_FM)*kron(idx_FM,speye(N_FM(2)))...
        +(t_FM(1)+t_int(1)+2*t_FM(2)+param.energylist_FM(kindex)-mu_FM)*kron(intx_FM,speye(N_FM(2)));
    
    
    
    
    K=[K_FM,K_int;K_int',K_Al];
    

    
%     K_bdg=[K,zero,zero,zero;
%            zero,K,zero,zero;
%            zero,zero,-K,zero;
%            zero,zero,zero,-K];

    K_bdg=kron(sparse([1],[1],[1],4,4),K)+...
        kron(sparse([2],[2],[1],4,4),K)+...
        kron(sparse([3],[3],[1],4,4),-K)+...
        kron(sparse([4],[4],[1],4,4),-K);
    
       
   
    H_bdg=K_bdg+D_bdg+Z_bdg;
    
   
    H_bdg=real((H_bdg+H_bdg')/2);
    k=max(20,2*ceil(1.05*k/2));
    if isempty(wfall{kindex})
        [vec,val]=eigs(H_bdg,k,'sm');
    else
        [vec,val]=eigs(H_bdg,k,'sm','StartVector',wfall{kindex}(:,1));
    end
        
    val=real(diag(val));
    if param.verbose==1
        fprintf('kindex=%d,min(val)=%f,max(val)=%f,',kindex,min(abs(val)),max(abs(val)));
    end

    
    %restart with larger k   
    while max(val)<param.ED
        if param.verbose==1
            fprintf('k (%d) is too small, restart with k (%d)\n',k,k*2);
        end
        k=k*2;
        [vec,val]=eigs(H_bdg,k,'sm');
        val=real(diag(val));
    end
    
    debyeindex=(abs(val)<=param.ED);
    val=val(debyeindex);
    vec=vec(:,debyeindex);
    [val,I]=sort(val);
    vec=vec(:,I);
    if param.verbose==1
        fprintf("len(val)=%d\n",length(val));
    end
    
    energyallknz{knzindex}=val;
    wfallknz{knzindex}=vec;

end

for knzindex=1:length(knzlist)
    
    energyall{knzlist(knzindex)}=energyallknz{knzindex};
    wfall{knzlist(knzindex)}=wfallknz{knzindex};
    energyall{end-knzlist(knzindex)+1}=energyall{knzlist(knzindex)};
    wfall{end-knzlist(knzindex)+1}=wfall{knzlist(knzindex)};
end
% toc;
end





