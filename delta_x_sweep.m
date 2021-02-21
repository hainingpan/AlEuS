xlist=linspace(0.2,0.23,60);
deltalist_s=xlist*0;
deltalist_n=xlist*0;
htotlist_s=xlist*0;
htotlist_n=xlist*0;

parfor xindex=1:length(xlist)
    [deltalist_s(xindex),htotlist_s(xindex)]=delta_x(xlist(xindex),0.2);
    [deltalist_n(xindex),htotlist_n(xindex)]=delta_x(xlist(xindex),0);
end
% figure;
% plot(xlist,deltalist);
% plot(1./xlist,log(deltalist),'-o')
% plot(xlist,deltalist_s,'DisplayName','$\Delta$');
% % legend;
figure;
hold on;
% plot(xlist,htotlist_s,'DisplayName','H_s');
% plot(xlist,htotlist_n,'DisplayName','H_n');
plot(xlist,htotlist_s-htotlist_n,'o-','DisplayName','H_s-H_n');
title(strcat('\Delta=',num2str(deltalist_s(1))))
% % legend;


function [d,htot]=delta_x(Ez,init)
param=main('Ez',Ez,'Nk',101,'mu',1,'T',0,'g',1);
kindex=1:param.Nk;
ux=(2*kindex-param.Nk-1)/(2*param.Nk);
klist=ux*param.b;
energylist=epsilon(klist,param);
param.energylist=energylist;

[energyall,wfall]=energyMF(init,param);

% figure;
dlist=[];
for i=1:1000
ave1=ave(energyall,wfall,param);
[energyall,wfall]=energyMF(ave1,param);
% plotband;
d=ave1*param.g;
% d=delta(energyall,wfall,param);
dlist=[dlist,d];
htot=totalenergy(energyall,wfall,param);
% plot(dlist);
% drawnow;
if length(dlist)>1
    if abs(dlist(end)-dlist(end-1))<1e-10
        break
    end
end
end
end
