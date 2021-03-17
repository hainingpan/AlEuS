Tlist=linspace(0,0.02,20);
gaplist_SC=Tlist*0;
gaplist_N=Tlist*0;

htotlist_N=Tlist*0;
htotlist_SC=Tlist*0;

for index=1:length(Tlist)
    [gaplist_N(index),htotlist_N(index)]=ns(Tlist(index));    
    [gaplist_SC(index),htotlist_SC(index)]=bs(Tlist(index),4e-3);

end
figure;
plot(Tlist,gaplist_SC-gaplist_N);


function [re,htot]=bs(Ez,gap)
param=main_3D('g',10,'N',[50,50,50],'ED',50*433*8.617333262e-5,'Ez',Ez,'pbc',1);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);

% upper=4e-3;
upper=gap;
lower=0e-3;
% figure;
middlelist=[];
wfall=cell(param.N(3),1);
while upper-lower>1e-8
    middle=(upper+lower)/2;
    ave0=middle*ones(param.N(1)*param.N(2),1);    
    [energyall,wfall]=energyMF_3D(ave0,wfall,param);
    ave1=ave_3D(energyall,wfall,param);
    fprintf('[%e,%e]\n',param.g*lower,param.g*upper);
    if mean(ave1,'all')-mean(ave0,'all')>0
        lower=middle;
    else
        upper=middle;
    end
    middlelist=[middlelist,middle];
%     plot(middlelist);
    drawnow;
end
% title(sprintf('%e',middle));
htot=totalenergy_3D(energyall,wfall,ave1,param);
re=middle*param.g;
end


function [re,htot]=ns(Ez)
param=main_3D('g',10,'N',[50,50,50],'ED',50*433*8.617333262e-5,'Ez',Ez,'pbc',1);

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
end