function param=main(varargin)
p=inputParser;
addParameter(p,'a',3.28e-10*5.076e6); % (eV^-1)
addParameter(p,'m',1);  % m_e, where m_e=0.511MeV
addParameter(p,'mu',0);
addParameter(p,'Ez',0);
addParameter(p,'g',1);
addParameter(p,'Nk',41);
addParameter(p,'T',0);


parse(p,varargin{:});
param=struct('a',p.Results.a,'m',p.Results.m*0.511e6,'mu',p.Results.mu,...
    'Ez',p.Results.Ez,'g',p.Results.g,'Nk',p.Results.Nk,'T',p.Results.T);
param.s0=[1,0;0,1];
param.sx=[0,1;1,0];
param.sy=[0,-1i;1i,0];
param.sz=[1,0;0,-1];
param.b=2*pi/param.a;
param.t=1/(2*param.m*param.a^2);


end


