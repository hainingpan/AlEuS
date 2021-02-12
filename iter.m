% param=main('Ez',0.5,'Nk',101,'mu',2,'T',0,'g',3,'a',3.28e-10*5.076e6);
param=main('Ez',0.1,'Nk',801,'mu',0.3,'T',0,'g',1,'a',3.28e-10*5.076e6);

kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(.0,param);

figure;
dlist=[];
hlist=[];
for i=1:1000
ave1=ave(energyall,wfall,param);
[energyall,wfall]=energyMF(ave1,param);
% plotband;
d=delta(energyall,wfall,param);
dlist=[dlist,d];
htot=totalenergy(energyall,wfall,param);
hlist=[hlist,htot];
% plot(dlist);
plot(hlist)
title(strcat('Delta=',num2str((d)),' E=',num2str(htot)));
drawnow;
if length(dlist)>1
    if abs(dlist(end)-dlist(end-1))<1e-10
        break
    end
end
end

