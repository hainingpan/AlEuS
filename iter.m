param=main('Ez',0.5,'Nk',101,'mu',2,'T',0,'g',3,'a',3.28e-10*5.076e6);
kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(0.2,param);

figure;
dlist=[];
for i=1:1000
ave1=ave(energyall,wfall,param);
[energyall,wfall]=energyMF(ave1,param);
% plotband;
d=delta(energyall,wfall,param);
dlist=[dlist,d];
plot(dlist);
title(strcat('Delta=',num2str((d))));
drawnow;
if length(dlist)>1
    if abs(dlist(end)-dlist(end-1))<1e-10
        break
    end
end
end

