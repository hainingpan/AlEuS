function d=delta(energyall,wfall,param)
%obsolete
%wfall(k,level,:)
Nk=size(wfall,1);
mat=zeros(4);
mat(4,1)=1;
mat(3,2)=-1;
fd=1./(exp((energyall)./param.T)+1);
d=0;
for kindex=1:Nk
    vec=squeeze(wfall(kindex,:,:)); 
    vec=vec.';
%     vec2=squeeze(wfall(kindex,2,:));
    d=d+param.g*(vec(:,1)'*mat*vec(:,1)*fd(kindex,1)...
                +vec(:,2)'*mat*vec(:,2)*fd(kindex,2)...
                +vec(:,3)'*mat*vec(:,3)*fd(kindex,3)...
                +vec(:,4)'*mat*vec(:,4)*fd(kindex,4));
%     d=d+param.g*(vec(:,1)'*mat*vec(:,1)+vec(:,2)'*mat*vec(:,2));
%      d=d+param.g*(conj(vec1(4))*vec1(1)+conj(vec2(4))*vec2(1));
end
d=d/(2*Nk);
% d=d/(Nk);
end
