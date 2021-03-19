function k=estimate_k_bilayer(kz,param)
% there is a problem here because the kx,ky is clearly not periodic
% Needs to find a better lower bound
kxindex=1:param.N_Al(1);
ux=(2*kxindex-param.N_Al(1)-1)/(2*param.N_Al(1));
kxlist=ux*param.b(1);

kyindex=1:param.N_Al(2);
uy=(2*kyindex-param.N_Al(2)-1)/(2*param.N_Al(2));
kylist=uy*param.b(2);

[kxgrid,kygrid]=meshgrid(kxlist,kylist);
Emesh=2*param.t_Al(1)*(1-cos(kxgrid*param.a(1)))+2*param.t_Al(2)*(1-cos(kygrid*param.a(2)))+2*param.t_Al(3)*(1-cos(kz*param.a(3)))-param.mu_Al;
k=4*nnz(Emesh<=param.ED & Emesh>=-param.ED)+0;
end