function R=thy(E)
mp=1;
mn=0.5;
kp=sqrt(2*mp*E);
kn=sqrt(2*mn*E);

B=(kn-kp)/(kn+kp);
% C=(2*kn)/(kn+kp);
R=abs(B)^2;
return

