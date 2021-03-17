function ham=Htb()
a=1;
mp=1;
mn=0.5;
Np=200;
Nn=200;
tp=1/(2*a^2*mp);
tn=1/(2*a^2*mn);

Mnn=2*tn*speye(Nn);
Mpp=2*tp*speye(Np);

Mnn=Mnn-tn*spdiags([ones(Nn,1),ones(Nn,1)],[-1,1],Nn,Nn);
Mpp=Mpp-tp*spdiags([ones(Np,1),ones(Np,1)],[-1,1],Np,Np);

Mnp=zeros(Nn,Np);
Mpn=zeros(Np,Nn);

Mnp(end,1)=3/2*tp-1/2*tn;
Mpn(1,end)=3/2*tn-1/2*tp;

ham=[Mnn,Mnp;Mpn,Mpp];
end
