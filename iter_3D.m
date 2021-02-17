param=main_3D('g',10,'N',[25,25,50]);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);


ave0=0.01/4*ones(param.N(1)*param.N(2),1);
[energyall,wfall]=energyMF_3D(ave0,param);

figure;
dlist=[];
for i=1:1000
    ave1=ave_3D(energyall,wfall,param);
    [energyall,wfall]=energyMF_3D(ave1,param);
    d=reshape(ave1*param.g,param.N(1),param.N(2));
    
    dlist(:,:,i)=d;
    subplot(3,1,1);
    plot(squeeze(mean(dlist,[1,2])));
    subplot(3,1,[2,3]);
    imagesc(d);
    daspect([1,1,1]);
    title(strcat('mean Delta=',num2str(mean(d,'all'))));
    colorbar;
    drawnow;
    if size(dlist,3)>1
        if abs(mean(dlist(:,:,end),'all')-mean(dlist(:,:,end-1),'all'))<1e-8
            break
        end
    end
end
    
    