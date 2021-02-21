function re=totalenergy_3D(energyall,wfall,ave,param)
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;
N=param.N;
t=param.t;
mu=param.mu;

delta=param.g*ave;
% Kinetic=2*sum(t)-mu+sum(param.muVar)/(N(1)*N(2));
Kinetic=0;
D=1/(N(1)*N(2))*sum(delta.^2)/param.g;

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
    

aveHbdg=0;
for i=1:size(energyall,1)
    if ~isempty(energyall{i})
        fd=1./(exp(energyall{i}./param.T)+1);
        wf=squeeze(wfall{i});
        K=-kron(Bx,speye(N(2)))*t(1)-kron(speye(N(1)),By)*t(2)+(2*t(1)+2*t(2)+param.energylist(i)-mu)*speye(N(1)*N(2))...
        +spdiags(param.muVar,0,N(1)*N(2),N(1)*N(2));
        K_bdg=[K,zero,zero,zero;
           zero,K,zero,zero;
           zero,zero,-K,zero;
           zero,zero,zero,-K];
       
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