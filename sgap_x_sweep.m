xlist=linspace(0,2,20);
deltalist=xlist*0;
parfor xindex=1:length(xlist)
    deltalist(xindex)=delta_x(xlist(xindex));
end
figure;
plot(xlist,deltalist);
% plot(1./xlist,log(deltalist),'-o')

function sgap=delta_x(Ez)
param=main('Ez',Ez,'Nk',401,'mu',2,'T',0,'g',3);
kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(1,param);

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

sgap=spectrum_gap(energyall);

end

