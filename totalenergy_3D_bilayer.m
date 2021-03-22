function re=totalenergy_3D_bilayer(energyall,wfall,ave,param)
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
N=param.N;
N_Al=param.N_Al;
N_FM=param.N_FM;
t_Al=param.t_Al;
t_FM=param.t_FM;
t_int=param.t_int;

mu_Al=param.mu_Al;
mu_FM=param.mu_FM;

delta=param.g*ave;
% Kinetic=2*sum(t)-mu+sum(param.muVar)/(N(1)*N(2));
Kinetic=0;
D=1/(N(1)*N(2))*sum(delta.^2)/param.g; %should I use N_Al or N


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
     
D_bdg=kron(sparse([1],[4],[1],4,4),-Dmat)...
    +kron(sparse([2],[3],[1],4,4),Dmat)...
    +kron(sparse([3],[2],[1],4,4),conj(Dmat))...
    +kron(sparse([4],[1],[1],4,4),-conj(Dmat));
     
Z_bdg=kron(kron(sz,sz),spdiags([ones(N_FM(1)*N_FM(2),1);zeros(N_Al(1)*N_Al(2),1)],0,N(1)*N(2),N(1)*N(2)))*param.Ez;
    

aveHbdg=0;
for i=1:size(energyall,1)
    if ~isempty(energyall{i})
        fd=1./(exp(energyall{i}./param.T)+1);
        wf=squeeze(wfall{i});
        
        K_Al=-kron(Bx_Al,speye(N_Al(2)))*t_Al(1)-kron(speye(N_Al(1)),By_Al)*t_Al(2)...
        +(2*t_Al(1)+2*t_Al(2)+param.energylist_Al(i)-mu_Al)*kron(idx_Al,speye(N_Al(2)))...
        +(t_Al(1)+t_int(1)+2*t_Al(2)+param.energylist_Al(i)-mu_Al)*kron(intx_Al,speye(N_Al(2)))...
        +spdiags(param.muVar,0,N_Al(1)*N_Al(2),N_Al(1)*N_Al(2));
    
        K_FM=-kron(Bx_FM,speye(N_FM(2)))*t_FM(1)-kron(speye(N_FM(1)),By_FM)*t_FM(2)...
            +(2*t_FM(1)+2*t_FM(2)+param.energylist_FM(i)-mu_FM)*kron(idx_FM,speye(N_FM(2)))...
            +(t_FM(1)+t_int(1)+2*t_FM(2)+param.energylist_FM(i)-mu_FM)*kron(intx_FM,speye(N_FM(2)));   
        
        K=[K_FM,K_int;K_int',K_Al];
              
        K_bdg=kron(sparse([1],[1],[1],4,4),K)+...
        kron(sparse([2],[2],[1],4,4),K)+...
        kron(sparse([3],[3],[1],4,4),-K)+...
        kron(sparse([4],[4],[1],4,4),-K);
    
        H_bdg=K_bdg+D_bdg+Z_bdg;
        
        H_bdg=real((H_bdg+H_bdg')/2);
       
        for j=1:size(energyall{i})
            aveHbdg=aveHbdg+wf(:,j)'*H_bdg*wf(:,j)*fd(j);
        end
    end
end
aveHbdg=aveHbdg/(2*prod(param.N));
        


re=aveHbdg+Kinetic+D;
end