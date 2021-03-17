param=main('Ez',0.1,'Nk',201,'mu',0.3,'T',0,'g',1,'a',3.28e-10*5.076e6,'m',1);
philist=linspace(-pi,pi,100);
Elist=1/(2*param.m*param.a^2)*philist.^2-param.mu;
figure;
subplot(1,2,1);
% hold on;
% plot(philist,Elist);
% plot(philist,-Elist);
subplot(1,2,2);
hold on;
Elist=2*(1/(2*param.m*param.a^2))*(1-cos(philist))-param.mu;
plot(philist,Elist+param.Ez);
plot(philist,Elist-param.Ez);

plot(philist,-(Elist+param.Ez));
plot(philist,-(Elist-param.Ez));


