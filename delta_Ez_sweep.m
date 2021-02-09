Ezlist=linspace(0,0.11,100);
deltalist=Ezlist*0;
parfor Ezindex=1:length(Ezlist)
    deltalist(Ezindex)=delta_Ez(Ezlist(Ezindex));
end
figure;
plot(Ezlist,deltalist);

function d=delta_Ez(Ez)
param=main('Ez',Ez,'Nk',101,'mu',0.5,'T',0,'g',0.5);
kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(1,param);

% figure;
dlist=[];
for i=1:1000
d=delta(energyall,wfall,param);
[energyall,wfall]=energyMF(d,param);
% plotband;
dlist=[dlist,d];

if length(dlist)>1
    if abs(dlist(end)-dlist(end-1))<1e-10
        break
    end
end
end
end

