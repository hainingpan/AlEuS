function param=main_3D_bilayer(varargin)
%For the bilayer
p=inputParser;
addParameter(p,'a',[1e-10*5.076e6,1e-10*5.076e6,1e-10*5.076e6]); % (eV^-1)
addParameter(p,'m_Al',1);  % m_e, where m_e=0.511MeV
addParameter(p,'m_FM',0.3);
addParameter(p,'mu_Al',11); % 11 eV in Al
addParameter(p,'mu_FM',-0.5); % 0.5 eV in FM

addParameter(p,'Ez',0); 
addParameter(p,'g',5.6);
% addParameter(p,'g_FM',0);
addParameter(p,'N_Al',[100,500,50]);
addParameter(p,'N_FM',[20,500,50]);
addParameter(p,'T',0);
addParameter(p,'ED',433*8.617333262e-5) % Debye frequency of Al
addParameter(p,'var',0.); %variance of disorder is 0.2eV
addParameter(p,'pbc',0); %variance of disorder is 0.2eV
addParameter(p,'verbose',1); % 0:suppress the output of each kz



parse(p,varargin{:});
param=struct('a',p.Results.a,'m_Al',p.Results.m_Al*0.511e6,'m_FM',p.Results.m_FM*0.511e6,...
    'mu_Al',p.Results.mu_Al,'mu_FM',p.Results.mu_FM,...
    'Ez',p.Results.Ez,'g',p.Results.g,...
    'N_Al',p.Results.N_Al,'N_FM',p.Results.N_FM,...
    'T',p.Results.T,'ED',p.Results.ED,'var',p.Results.var,'pbc',p.Results.pbc,...
    'verbose',p.Results.verbose);
param.s0=[1,0;0,1];
param.sx=[0,1;1,0];
param.sy=[0,-1i;1i,0];
param.sz=[1,0;0,-1];
param.b=2*pi./param.a; 
param.t_Al=1./(2*param.m_Al*param.a.^2);
param.t_FM=1./(2*param.m_FM*param.a.^2);
param.t_int=1./(2*sqrt(param.m_Al*param.m_FM)*param.a.^2);
param.N=param.N_Al;
param.N(1)=param.N_Al(1)+param.N_FM(1);
% param.g=[param.g_FM*ones(param.N_FM(1)*param.N_FM(2)),param.g_Al*ones(param.N_Al(1)*param.N_Al(2)];
% param.muVar=param.var*randn(param.N_Al(1)*param.N_Al(2),1);
param.muVar=param.var*(2*rand(param.N_Al(1)*param.N_Al(2),1)-1);



end


