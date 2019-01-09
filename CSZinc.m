%% Cargar datos
% Secciones de EXFOR.
load('CrossSections.mat');
PG_Zn64_L1_CS=PG_Zn64_L1_CS./1000;
PG_Zn64_L2_CS=PG_Zn64_L2_CS./1000;
PG_Zn64_L3_CS=PG_Zn64_L3_CS./1000;
PG_Zn64_L4_CS=PG_Zn64_L4_CS./1000;
PG_Zn64_L5_CS=PG_Zn64_L5_CS./1000;
PG_Zn66_L1_CS=PG_Zn66_L1_CS./1000;
PG_Zn66_L2_CS=PG_Zn66_L2_CS./1000;
PG_Zn66_L3_CS=PG_Zn66_L3_CS./1000;
PG_Zn66_L4_CS=PG_Zn66_L4_CS./1000;
PG_Zn68_L1_CS=PG_Zn68_L1_CS./1000;
PG_Zn68_L2_CS=PG_Zn68_L2_CS./1000;
PG_Zn68_L3_CS=PG_Zn68_L3_CS./1000;
PG_Zn68_L4_CS=PG_Zn68_L4_CS./1000;
PG_Zn68_L5_CS=PG_Zn68_L5_CS./1000;

Zn64_ab = 0.492; %4
Zn66_ab = 0.277; %5
Zn68_ab = 0.185; %6
Ca44_ab = 0.0223233841; %

%% Fit L1 Zn 64

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn64_L1_CS,'bo')
hold on
Zn64_L1_F = fit(PG_Zn64_E,PG_Zn64_L1_CS,'smoothingspline','SmoothingParam',0.99);
Zn64_L1_F.p.coefs(1,:) = [0 0 0 0];
Zn64_L1_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn64_L1_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L2 Zn 64

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn64_L2_CS,'bo')
hold on
Zn64_L2_F = fit(PG_Zn64_E,PG_Zn64_L2_CS,'smoothingspline','SmoothingParam',0.99);
Zn64_L2_F.p.coefs(1,:) = [0 0 0 0];
Zn64_L2_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn64_L2_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L3 Zn 64

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn64_L3_CS,'bo')
hold on
Zn64_L3_F = fit(PG_Zn64_E,PG_Zn64_L3_CS,'smoothingspline','SmoothingParam',0.99);
Zn64_L3_F.p.coefs(1,:) = [0 0 0 0];
Zn64_L3_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn64_L3_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L4 Zn 64

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn64_L4_CS,'bo')
hold on
Zn64_L4_F = fit(PG_Zn64_E,PG_Zn64_L4_CS,'smoothingspline','SmoothingParam',0.99);
Zn64_L4_F.p.coefs(1,:) = [0 0 0 0];
Zn64_L4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn64_L4_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L5 Zn 64

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn64_L5_CS,'bo')
hold on
Zn64_L5_F = fit(PG_Zn64_E,PG_Zn64_L5_CS,'smoothingspline','SmoothingParam',0.99);
Zn64_L5_F.p.coefs(1,:) = [0 0 0 0];
Zn64_L5_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn64_L5_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L1 Zn 66

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn66_L1_CS,'bo')
hold on
Zn66_L1_F = fit(PG_Zn64_E,PG_Zn66_L1_CS,'smoothingspline','SmoothingParam',0.99);
Zn66_L1_F.p.coefs(1,:) = [0 0 0 0];
Zn66_L1_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn66_L1_F(Eval),'b-');
axis([0 100 0 0.2]);



%% Fit L2 Zn 66

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn66_L2_CS,'bo')
hold on
Zn66_L2_F = fit(PG_Zn64_E,PG_Zn66_L2_CS,'smoothingspline','SmoothingParam',0.99);
Zn66_L2_F.p.coefs(1,:) = [0 0 0 0];
Zn66_L2_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn66_L2_F(Eval),'b-');
axis([0 100 0 0.05]);

%% Fit L3 Zn 66

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn66_L3_CS,'bo')
hold on
Zn66_L3_F = fit(PG_Zn64_E,PG_Zn66_L3_CS,'smoothingspline','SmoothingParam',0.99);
Zn66_L3_F.p.coefs(1,:) = [0 0 0 0];
Zn66_L3_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn66_L3_F(Eval),'b-');
axis([0 100 0 0.2]);


%% Fit L4 Zn 66

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn66_L4_CS,'bo')
hold on
Zn66_L4_F = fit(PG_Zn64_E,PG_Zn66_L4_CS,'smoothingspline','SmoothingParam',0.99);
Zn66_L4_F.p.coefs(1,:) = [0 0 0 0];
Zn66_L4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn66_L4_F(Eval),'b-');
axis([0 100 0 0.2]);

%% Fit L1 Zn 68

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn68_L1_CS,'bo')
hold on
Zn68_L1_F = fit(PG_Zn64_E,PG_Zn68_L1_CS,'smoothingspline','SmoothingParam',0.99);
Zn68_L1_F.p.coefs(1,:) = [0 0 0 0];
Zn68_L1_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn68_L1_F(Eval),'b-');
axis([0 100 0 0.05]);


%% Fit L2 Zn 68

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn68_L2_CS,'bo')
hold on
Zn68_L2_F = fit(PG_Zn64_E,PG_Zn68_L2_CS,'smoothingspline','SmoothingParam',0.99);
Zn68_L2_F.p.coefs(1,:) = [0 0 0 0];
Zn68_L2_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn68_L2_F(Eval),'b-');
axis([0 100 0 0.05]);



%% Fit L3 Zn 68

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn68_L3_CS,'bo')
hold on
Zn68_L3_F = fit(PG_Zn64_E,PG_Zn68_L3_CS,'smoothingspline','SmoothingParam',0.99);
Zn68_L3_F.p.coefs(1,:) = [0 0 0 0];
Zn68_L3_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn68_L3_F(Eval),'b-');
axis([0 100 0 0.05]);


%% Fit L4 Zn 68

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn68_L4_CS,'bo')
hold on
Zn68_L4_F = fit(PG_Zn64_E,PG_Zn68_L4_CS,'smoothingspline','SmoothingParam',0.99);
Zn68_L4_F.p.coefs(1,:) = [0 0 0 0];
Zn68_L4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn68_L4_F(Eval),'b-');
axis([0 100 0 0.05]);


%% Fit L5 Zn 68

figure
Eval = 0:0.1:300; % MeV
plot(PG_Zn64_E,PG_Zn68_L5_CS,'bo')
hold on
Zn68_L5_F = fit(PG_Zn64_E,PG_Zn68_L5_CS,'smoothingspline','SmoothingParam',0.99);
Zn68_L5_F.p.coefs(1,:) = [0 0 0 0];
Zn68_L5_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Zn68_L5_F(Eval),'b-');
axis([0 100 0 0.05]);
