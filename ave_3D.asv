function ave=ave_3D(energyall,wfall,param)
% 1/2 1/Nz*sum_kz c_{r,-k,down} c_{r,k,up}
N=param.N;
% mat=zeros(4);
% mat(4,1)=1;
% mat(3,2)=-1;
% debye=(abs(energyall)<=param.ED);
ave=zeros(N(1)*N(2),1);
for i=1:size(energyall,1)
    if ~isempty(energyall{i})
        fd=1./(exp(energyall{i}./param.T)+1);
        wf=squeeze(wfall{i}); % 4N*4N
        wf3=reshape(wf,4,N(1)*N(2),size(wf,2));
        ave=ave+squeeze(conj(wf3(4,:,:)).*wf3(1,:,:)-conj(wf3(3,:,:)).*wf3(2,:,:))*fd;
    end
end
ave=ave/(2*N(3));
end