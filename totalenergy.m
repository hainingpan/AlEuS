function re=totalenergy(energyall,wfall,param)
sz=param.sz;
sy=param.sy;
sx=param.sx;
s0=param.s0;

d=delta(energyall,wfall,param);

% fd_K=1./(exp((param.energylist)./param.T)+1);
% K=1/param.Nk*fd_K*param.energylist';
K=1/param.Nk*sum(param.energylist');


fd=1./(exp((energyall)./param.T)+1);
aveHbdg=0;
for kindex=1:param.Nk

    Hbdg=param.energylist(kindex)*kron(sz,s0)+d*kron(sy,sy)+param.Ez*kron(sz,sz);
    vec=squeeze(wfall(kindex,:,:)); 
    vec=vec.';
    aveHbdg=aveHbdg+(vec(:,1)'*Hbdg*vec(:,1)*fd(kindex,1)...
                +vec(:,2)'*Hbdg*vec(:,2)*fd(kindex,2)...
                +vec(:,3)'*Hbdg*vec(:,3)*fd(kindex,3)...
                +vec(:,4)'*Hbdg*vec(:,4)*fd(kindex,4));
end
aveHbdg=aveHbdg/param.Nk/2;
re=d^2/param.g+K+aveHbdg;
end
