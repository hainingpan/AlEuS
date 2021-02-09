param=main('Ez',0.00,'Nk',100,'mu',0.5,'T',0.0,'g',0.3);
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

