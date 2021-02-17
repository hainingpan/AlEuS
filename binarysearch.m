Tlist=linspace(0,0.02,50);
gaplist=Tlist*0;
for index=1:length(Tlist)
    gaplist(index)=bs(Tlist(index));
end
figure;
plot(Tlist,gaplist);


function re=bs(T)
param=main_3D('g',10,'N',[25,25,50],'ED',50*433*8.617333262e-5,'T',T);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3);
param.energylist=epsilon_3D(param.klist,param);

upper=4e-3;
lower=0e-3;
% figure;
middlelist=[];
while upper-lower>1e-6
    middle=(upper+lower)/2;
    ave0=middle*ones(param.N(1)*param.N(2),1);    
    [energyall,wfall]=energyMF_3D(ave0,param);
    ave1=ave_3D(energyall,wfall,param);
    fprintf('[%e,%e]\n',lower,upper);
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
re=middle*param.g;
end