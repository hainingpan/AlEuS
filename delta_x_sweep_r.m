xlist=linspace(0.3,0.9,20);
deltalist=xlist*0;
parfor xindex=1:length(xlist)
    deltalist(xindex)=delta_x(xlist(xindex));
end
figure;
% plot(xlist,deltalist);
plot(1./xlist,log(deltalist),'-o')

function dm=delta_x(g)
param=main('Ez',0.,'Nk',100,'mu',0.5,'T',0,'g',g);
ave0=50*ones(param.Nk,1);
[energyall,wfall]=energyMF_r(ave0,param);

dlist=[];
for i=1:1000
    ave1=ave_r(energyall,wfall,param);
    [energyall,wfall]=energyMF_r(ave1,param);
    % plotband;
    % d=delta(energyall,wfall,param);
    d=ave1*param.g/2;
    dlist=[dlist,d];
    title(strcat('mean Delta=',num2str(mean(d))));
    drawnow;
    if size(dlist,2)>1
        if abs(mean(dlist(:,end))-mean(dlist(:,end-1)))<1e-5
            break
        end
    end
end
dm=mean(d);
end