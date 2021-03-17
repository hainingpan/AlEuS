ham=Htb();
[V,D]=eig(full(ham));
D=diag(D);
[~,I]=sort(D);
V=V(:,I);

figure;
level=1;
plot(abs(V(:,level)).^2);
sum(V(1:end/2,level).^2)



