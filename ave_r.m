function ave=ave_r(energyall,wfall,param)
%<c_{i,up}c_{i,down}>
Nk=param.Nk;

mat=zeros(4);
mat(4,1)=-1;
mat(3,2)=1;
fd=1./(exp((energyall)./param.T)+1);
ave=zeros(Nk,1);

vec=wfall;
for i=1:Nk    
    for j=1:4*Nk
        ave(i)=ave(i)+vec([i,Nk+i,2*Nk+i,3*Nk+i],j)'*mat*vec([i,Nk+i,2*Nk+i,3*Nk+i],j)*fd(j);
    end
end

end
