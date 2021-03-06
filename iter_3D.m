param=main_3D('g',5.55,'N',[60,100,50],'ED',50*433*8.617333262e-5,'pbc',0,'Ez',0,'verbose',0);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3); 
param.energylist=epsilon_3D(param.klist,param);


% ave0=0.0016*ones(param.N(1)*param.N(2),1);
ave0=2e-3*ones(param.N(1)*param.N(2),1);

[energyall,wfall]=energyMF_3D(ave0,cell(param.N(3),1),param);

figure;
dlist=[];
htotlist=[];
for i=1:1000
    ave1=ave_3D(energyall,wfall,param);
    d=reshape(ave1*param.g,param.N(2),param.N(1));
    dlist(:,:,i)=d;
    htot=totalenergy_3D(energyall,wfall,ave1,param);
    htotlist=[htotlist,htot];
    subplot(3,2,1);
    plot(squeeze(mean(dlist,[1,2])));
    ylabel('mean Delta')
    subplot(3,2,2);
    plot(htotlist);
    ylabel('total energy')
    subplot(3,2,[3,4,5,6]);
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
    [energyall,wfall]=energyMF_3D(ave1,wfall,param);
end
    
    