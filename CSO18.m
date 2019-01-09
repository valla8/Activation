%% Cargar datos
% Secciones de EXFOR.

load('CrossSections.mat');
PG_O18_L1_CS=PG_O18_L1_CS./1000;
PG_O18_L2_CS=PG_O18_L2_CS./1000;
PG_O18_L3_CS=PG_O18_L3_CS./1000;
PG_O18_L4_CS=PG_O18_L4_CS./1000;
PG_O18_L5_CS=PG_O18_L5_CS./1000;


%% Fit O18 L1
figure
Eval = 0:0.1:300; % MeV
plot(PG_O18_E,PG_O18_L1_CS,'bo')
hold on
O18_L1_F = fit(PG_O18_E,PG_O18_L1_CS,'smoothingspline','SmoothingParam',0.99);
O18_L1_F.p.coefs(1,:) = [0 0 0 0];
O18_L1_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_L1_F(Eval),'b-');
axis([0 100 0 0.2]);



%% Fit O18 L2
figure
Eval = 0:0.1:300; % MeV
plot(PG_O18_E,PG_O18_L2_CS,'bo')
hold on
O18_L2_F = fit(PG_O18_E,PG_O18_L2_CS,'smoothingspline','SmoothingParam',0.99);
O18_L2_F.p.coefs(1,:) = [0 0 0 0];
O18_L2_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_L2_F(Eval),'b-');
axis([0 100 0 0.2]);





%% Fit O18 L3
figure
Eval = 0:0.1:300; % MeV
plot(PG_O18_E,PG_O18_L3_CS,'bo')
hold on
O18_L3_F = fit(PG_O18_E,PG_O18_L3_CS,'smoothingspline','SmoothingParam',0.99);
O18_L3_F.p.coefs(1,:) = [0 0 0 0];
O18_L3_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_L3_F(Eval),'b-');
axis([0 100 0 0.2]);


%% Fit O18 L4
figure
Eval = 0:0.1:300; % MeV
plot(PG_O18_E,PG_O18_L4_CS,'bo')
hold on
O18_L4_F = fit(PG_O18_E,PG_O18_L4_CS,'smoothingspline','SmoothingParam',0.99);
O18_L4_F.p.coefs(1,:) = [0 0 0 0];
O18_L4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_L4_F(Eval),'b-');
axis([0 100 0 0.2]);




%% Fit O18 L5
figure
Eval = 0:0.1:300; % MeV
plot(PG_O18_E,PG_O18_L5_CS,'bo')
hold on
O18_L5_F = fit(PG_O18_E,PG_O18_L5_CS,'smoothingspline','SmoothingParam',0.99);
O18_L55F.p.coefs(1,:) = [0 0 0 0];
O18_L4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_L5_F(Eval),'b-');
axis([0 100 0 0.2]);


