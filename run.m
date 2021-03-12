function run(Nkz,g)
param=main_3D('g',g,'N',[60,100,Nkz],'ED',50*433*8.617333262e-5,'pbc',0,'Ez',0.0,'verbose',0);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);

ave0=2e-3*ones(param.N(1)*param.N(2),1);

[energyall,wfall]=energyMF_3D(ave0,cell(param.N(3),1),param);

dlist=[];
htotlist=[];
for i=1:1000
    ave1=ave_3D(energyall,wfall,param);
    d=reshape(ave1*param.g,param.N(1),param.N(2));
    dlist(:,:,i)=d;
    htot=totalenergy_3D(energyall,wfall,ave1,param);
    htotlist=[htotlist,htot];   
    fprintf('Average Gap: %e (meV), Total energy: %e (meV)\n',mean(d),mean(htot));
    if size(dlist,3)>1
        if abs(mean(dlist(:,:,end),'all')-mean(dlist(:,:,end-1),'all'))<1e-8
            break
        end
    end
    
    save(sprintf('Nk(%d,%d,%d)g%dED%d.mat',param.N,param.g,param.ED/(433*8.617333262e-5)),'dlist','htotlist','param');
    [energyall,wfall]=energyMF_3D(ave1,wfall,param);
end
end
    