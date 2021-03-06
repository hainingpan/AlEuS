function [Elist,dos]=DOS(energyall,param)
energyall_sort=sort(energyall(:));

eta=1e-4;
Emin=energyall_sort(1);
Emax=energyall_sort(end);
Elist=linspace(Emin,Emax,1000);
energyproj=energyall(:);
deltaf=1/pi*eta./((Elist-energyproj).^2+eta^2); %axis 1: state; axis 2: Elist
dos=sum(deltaf,1);

dos=dos/size(energyall,1)/(param.a/5.076e-3);  % of eV^-1 nm^-1
figure;plot(Elist,dos);
xlabel('E (meV)');
ylabel('DOS(eV^{-1}nm^{-1})');
end
