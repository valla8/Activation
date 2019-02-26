load('Efi-100m.mat')

plot(Energ,Efi,'bo');hold on;
FF=fit(Energ,Efi,'smoothingspline','SmoothingParam',0.2)
FF.p.coefs(1,:)=[0 0 0 0];
FF.p.coefs(end,:) = [0 0 0 0];
Eval=linspace(1,1000,1000);
plot(Eval,FF(Eval))