Tlist=linspace(0,0.1,100);
deltalist=Tlist*0;
parfor Tindex=1:length(Tlist)
    deltalist(Tindex)=delta_T(Tlist(Tindex));
end
figure;
plot(Tlist,deltalist);

function d=delta_T(T)
param=main('Ez',0,'Nk',101,'mu',0.5,'T',T,'g',0.5);
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

