function [energy_Al,energy_FM]=epsilon_3D_bilayer(k,param)
energy_Al=2*param.t_Al(3)*(1-cos(param.a(3)*k));
energy_FM=2*param.t_FM(3)*(1-cos(param.a(3)*k));

end