Tlist=linspace(0,0.15,50);
deltalist=Tlist*0;
for Tindex=1:length(Tlist)
    deltalist(Tindex)=delta_T(Tlist(Tindex));
end
figure;
plot(Tlist,deltalist);

function d=delta_T(T)
param=main('Ez',0,'Nk',401,'mu',0.3,'T',T,'g',1,'ED',0.7);
kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(0.1,param);

% figure;
dlist=[];
for i=1:1000
ave1=ave(energyall,wfall,param);
[energyall,wfall]=energyMF(ave1,param);
% plotband;
d=ave1*param.g;
% d=delta(energyall,wfall,param);
dlist=[dlist,d];
% plot(dlist);
% title(strcat('Delta=',num2str((d))));
% drawnow;
if length(dlist)>1
    if abs(dlist(end)-dlist(end-1))<1e-10
        break
    end
end
end
end

