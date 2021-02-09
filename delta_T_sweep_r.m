Tlist=linspace(0,0.03,20);
deltalist=Tlist*0;
parfor Tindex=1:length(Tlist)
    deltalist(Tindex)=delta_T(Tlist(Tindex));
end
figure;
plot(Tlist,deltalist);

function dm=delta_T(T)
param=main('Ez',0.,'Nk',100,'mu',0.5,'T',T,'g',0.5);
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