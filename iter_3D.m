param=main_3D('g',1,'N',[100,100,800]);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);


ave0=0*ones(param.N(1)*param.N(2),1);
[energyall,wfall]=energyMF_3D(ave0,param);

figure;
dlist={};
for i=1:1000
    ave1=ave_3D(energyall,wfall,param);
    [energyall,wfall]=energyMF_3D(ave1,param);
    d=reshape(ave1*param.g,param.N(1),param.N(2));
    
    dlist{i}=d;
    imagesc(d);
    title(strcat('mean Delta=',num2str(mean(d,'all'))));
    colorbar;
    drawnow;
    if size(dlist,2)>1
        if abs(mean(dlist{end},'all')-mean(dlist{end-1},'all'))<1e-10
            break
        end
    end
end
    
    