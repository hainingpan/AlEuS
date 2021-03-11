%Verify the Vz dependence

Tlist=linspace(0,0.02,20);
gaplist_SC=Tlist*0;
gaplist_N=Tlist*0;

htotlist_N=Tlist*0;
htotlist_SC=Tlist*0;

for index=1:length(Tlist)
    fprintf("index=%d\n",index);
    [gaplist_N(index),htotlist_N(index)]=ns(Tlist(index));    
    [gaplist_SC(index),htotlist_SC(index)]=sc(Tlist(index));

end
figure;
plot(Tlist,gaplist_SC-gaplist_N);
save(sprintf("Zeeman.mat"),'Tlist','gaplist_SC','gaplist_N','htotlist_N','htotlist_SC');

function [re,htot]=sc(Ez)
param=main_3D('g',10,'N',[50,50,50],'ED',50*433*8.617333262e-5,'Ez',Ez,'pbc',0,'verbose',0);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);

ave0=1.87e-3*ones(param.N(1)*param.N(2),1);

[energyall,wfall]=energyMF_3D(ave0,cell(param.N(3),1),param);

dlist=[];
htotlist=[];
for i=1:1000
    ave1=ave_3D(energyall,wfall,param);
    d=reshape(ave1*param.g,param.N(1),param.N(2));
    dlist(:,:,i)=d;
    htot=totalenergy_3D(energyall,wfall,ave1,param);
    htotlist=[htotlist,htot];   
    fprintf('(SC) Average Gap: %e (meV), Total energy: %e (meV)\n',mean(d,'all'),mean(htot));
    if size(dlist,3)>1
        if abs(mean(dlist(:,:,end),'all')-mean(dlist(:,:,end-1),'all'))<1e-8
            break
        end
    end    
    [energyall,wfall]=energyMF_3D(ave1,wfall,param);
end
re=mean(d,'all');
end


function [re,htot]=ns(Ez)
param=main_3D('g',10,'N',[50,50,50],'ED',50*433*8.617333262e-5,'Ez',Ez,'pbc',0,'verbose',0);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);


wfall=cell(param.N(3),1);

ave0=0*ones(param.N(1)*param.N(2),1);    
[energyall,wfall]=energyMF_3D(ave0,wfall,param);
ave1=ave_3D(energyall,wfall,param);

htot=totalenergy_3D(energyall,wfall,ave1,param);
re=mean(ave1,'all')*param.g;

fprintf('(NS) Average Gap: %e (meV), Total energy: %e (meV)\n',re,mean(htot));

end