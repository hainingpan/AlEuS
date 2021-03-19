param=main_3D_bilayer('N_Al',[60,100,50],'N_FM',[1,100,50],'g',5.55,'ED',50*433*8.617333262e-5);

kindex=1:param.N(3);
uz=(2*kindex-param.N(3)-1)/(2*param.N(3));
param.klist=uz*param.b(3); 
param.energylist=epsilon_3D_bilayer(param.klist,param);

ave0=3.5e-4/param.g*[zeros(param.N_FM(1)*param.N_FM(2),1);ones(param.N_Al(1)*param.N_Al(2),1)];

[energyall,wfall]=energyMF_3D_bilayer(ave0,cell(param.N(3),1),param);
 
figure;
dlist=[];
htotlist=[];
for i=1:1000
    ave1=ave_3D_bilayer(energyall,wfall,param);
    d=reshape(ave1*param.g,param.N(2),param.N(1));
    dlist(:,:,i)=d;
%     htot=totalenergy_3D(energyall,wfall,ave1,param);
%     htotlist=[htotlist,htot];
    subplot(3,2,1);
    plot(squeeze(mean(dlist(:,param.N_FM(1)+1:end,:),[1,2])));
    ylabel('mean Delta')
    subplot(3,2,2);
    plot(htotlist);
    ylabel('total energy')
    subplot(3,2,[3,4,5,6]);
    imagesc(d);
    daspect([1,1,1]);
    title(strcat('mean Delta in Al=',num2str(mean(d(:,param.N_FM(1)+1:end),'all'))));
    colorbar;
    drawnow;
    if size(dlist,3)>1
        if abs(mean(dlist(:,param.N_FM(1)+1:end,end),'all')-mean(dlist(:,param.N_FM(1)+1:end,end-1),'all'))<1e-8
            break
        end
    end
    [energyall,wfall]=energyMF_3D_bilayer(ave1,wfall,param);
end
    
    

