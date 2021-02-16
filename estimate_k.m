function k=estimate_k(kz,param)
kxindex=1:param.N(1);
ux=(2*kxindex-param.N(1)-1)/(2*param.N(1));
kxlist=ux*param.b(1);

kyindex=1:param.N(2);
uy=(2*kyindex-param.N(2)-1)/(2*param.N(2));
kylist=uy*param.b(2);

[kxgrid,kygrid]=meshgrid(kxlist,kylist);
Emesh=2*param.t(1)*(1-cos(kxgrid*param.a(1)))+2*param.t(2)*(1-cos(kygrid*param.a(2)))+2*param.t(3)*(1-cos(kz*param.a(3)))-param.mu;
k=4*nnz(Emesh<=param.ED & Emesh>=-param.ED)+0;
end