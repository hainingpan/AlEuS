% param=main('Ez',0.00,'Nk',200,'mu',1,'T',0.1,'g',1,'a',3.28e-10*5.076e6);
param=main('Ez',0.0,'Nk',101,'mu',2,'T',0,'g',3,'a',3.28e-10*5.076e6);

ave0=50*ones(param.Nk,1);
[energyall,wfall]=energyMF_r(ave0,param);

figure;
dlist=[];
for i=1:1000
ave1=ave_r(energyall,wfall,param);
[energyall,wfall]=energyMF_r(ave1,param);
% plotband;
% d=delta(energyall,wfall,param);
d=ave1*param.g/2;
dlist=[dlist,d];
plot(d);
title(strcat('mean Delta=',num2str(mean(d))));
drawnow;
if size(dlist,2)>1
    if abs(mean(dlist(:,end))-mean(dlist(:,end-1)))<1e-5
        break
    end
end
end

