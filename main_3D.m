function param=main_3D(varargin)
p=inputParser;
addParameter(p,'a',[1e-10*5.076e6,1e-10*5.076e6,1e-10*5.076e6]); % (eV^-1)
addParameter(p,'m',1);  % m_e, where m_e=0.511MeV
addParameter(p,'mu',11); % 11 eV in Al
addParameter(p,'Ez',0); 
addParameter(p,'g',1);
addParameter(p,'N',[100,100,100]);
addParameter(p,'T',0);
addParameter(p,'ED',433*8.617333262e-5) % Debye frequency of Al
addParameter(p,'var',0.); %variance of disorder is 0.2eV
addParameter(p,'pbc',0); %variance of disorder is 0.2eV
addParameter(p,'verbose',1); % 0:suppress the output of each kz



parse(p,varargin{:});
param=struct('a',p.Results.a,'m',p.Results.m*0.511e6,'mu',p.Results.mu,...
    'Ez',p.Results.Ez,'g',p.Results.g,'N',p.Results.N,'T',p.Results.T,...
    'ED',p.Results.ED,'var',p.Results.var,'pbc',p.Results.pbc,'verbose',p.Results.verbose);
param.s0=[1,0;0,1];
param.sx=[0,1;1,0];
param.sy=[0,-1i;1i,0];
param.sz=[1,0;0,-1];
param.b=2*pi./param.a; 
param.t=1./(2*param.m*param.a.^2);
param.muVar=param.var*randn(param.N(1)*param.N(2),1);



end


