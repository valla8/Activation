% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es

%% Par?metros a modificar:
clear all
Zn_fraction = 0.1; % Fraccion de Zn (tanto por uno)
dx = 0.1; % Espaciado de la malla (cm)
E0 = 58; % Energ?a inicial del haz (MeV)

%% Cargar datos
% Secciones de EXFOR.
load('CrossSections2.mat');
%C12_C10_CS=C12_C10_CS./1000;
%C12_C11_CS=C12_C11_CS./1000;
%PG_C12_C12_4_CS=PG_C12_C12_4_CS./1000;
%PG_O16_C12_4_CS=PG_O16_C12_4_CS./1000;
%PG_N14_N14_1_CS=PG_N14_N14_1_CS./1000;
%PG_N14_N14_2_CS=PG_N14_N14_2_CS./1000;
%PG_O16_O16_6_CS=PG_O16_O16_6_CS./1000;
%PG_C12_C12_4_E=PG_C12_C12_4_E;
%PG_O16_C12_4_E=PG_O16_C12_4_E;
%PG_N14_N14_1_E=PG_N14_N14_1_E;
%PG_N14_N14_2_E=PG_N14_N14_2_E;
%PG_O16_O16_6_E=PG_O16_O16_6_E;

%PG Nitrogeno Tallys
%Para poner experimentales comentar estas 4 lineas
PG_N14_N14_2_CS=PG_N14_L2_CS./1000;
PG_N14_N14_1_CS=PG_N14_L1_CS./1000;
PG_N14_N14_1_E=PG_N14_L1_E;
PG_N14_N14_2_E=PG_N14_L2_E;

I127_Xe127_CS=I127_Xe127_CS./1000;
I127_Xe125_CS=I127_Xe125_CS./1000;
I127_Xe123_CS=I127_Xe123_CS./1000;
I127_Xe122_CS=I127_Xe122_CS./1000;
I127_Xem_CS=I127_Xem_CS./1000;

Na23_Mg23_CS=Na23_Mg23_CS./1000;



% Vidas medias en s
load('MeanLives.mat');
landa_C10 = log(2) / T_C10;
landa_C11 = log(2) / T_C11;
landa_Ga64 = log(2) / T_Ga64;
landa_Ga65 = log(2) / T_Ga65;
landa_Ga66 = log(2) / T_Ga66;
landa_Ga67 = log(2) / T_Ga67;
landa_Ga68 = log(2) / T_Ga68;
landa_N13 = log(2) / T_N13;
landa_O15 = log(2) / T_O15;
landa_Sc44 = log(2) / T_Sc44;
landa_Xe127 = log(2) / T_Xe127;
landa_Xe125 = log(2) / T_Xe125;
landa_Xe123 = log(2) / T_Xe123;
landa_Xe122 = log(2) / T_Xe122;
landa_Xe127m = log(2) / T_Xe127m;
landa_Mg23= log(2) / T_Mg23

% Natural abundances
O16_ab = 0.99729; %1
N14_ab = 0.996; %2
C12_ab = 0.989; %3
Zn64_ab = 0.492; %4
Zn66_ab = 0.277; %5
Zn68_ab = 0.185; %6
Zn67_ab = 0.04;
Ca44_ab = 0.0223233841; %
I127_ab = 1;
ab = [O16_ab N14_ab C12_ab Zn64_ab Zn66_ab Zn68_ab];

% Tissue compositions (fuente: Zhu 2011, ICRU Report #63, ICRU Report #44)
% Material H (%) C (%) N (%) O (%) Zn64(%) Zn66(%) Zn68(%)  Density (g ml?1)
% Soft tissue 0.102 0.143 0.034 0.708 0 0 0 1.06
% Water: 0.0305 0 0 0.9695 0 0 0 1.00
Comp_tissue = [0.619977 0.118053 0.0205881 0.240345 0 0 ];
Comp_tissue_2 = [0.619977 0.118053 0.0205881 0.240345 0 0 0];
Comp_water = [0.6687 0 0 0.3313 0 0 0];
Comp_Zn = [0 0 0 0 Zn64_ab Zn66_ab Zn68_ab];
Comp_tissue_Zn =Zn_fraction*Comp_Zn + (1-Zn_fraction)* Comp_tissue_2;!
Comp_water_Zn = Zn_fraction*Comp_Zn + (1-Zn_fraction)*Comp_water;
% Material H (%) C (%) N (%) O (%) P (%) Ca(%)  Density (g ml?1)
Comp_bone=[0.475402    0.121899    0.0304112    0.282845    0.0343791    0.0531337];
W_ele = [1.00794 12.011 14.00674 15.9994 30.973762 40.078];
Comp_adipose = [0.634657 0.284  0.0030464 0.0777462 0 0];
Comp_PMMA= [0.535002   0.332244 0  0.132754 0 0];
Comp_PMMA_w= [0.080   0.60 0  0.32 0 0];

% Cargar stopping powers (only for water, tissue, Zn, bone)
load('stoppingpowers.mat');

%% Fit secciones eficaces ZN-Ga
figure
Eval = 0:0.1:300; % MeV
hold off
plot(Zn66_Ga66_E,Zn66_Ga66_CS,'bo')
hold on
Zn66_Ga66_F = fit(Zn66_Ga66_E,Zn66_Ga66_CS,'smoothingspline','SmoothingParam',0.1);
Zn66_Ga66_F.p.coefs(1,:) = [0 0 0 0];
Zn66_Ga66_F.p.coefs(end,:) = [0 0 0 0];

axis([0 100 0 0.8]);

% Las secciones eficaces de Zn64 y Zn68 se calculan modificando
% manualmente el fit para la secci?n eficaz del Zn66 ajustando
% a los puntos experimentales existentes.

plot(Eval,Zn66_Ga66_F(Eval),'b-');
plot(Zn64_Ga64_E,Zn64_Ga64_CS,'xg')
Zn64_Ga64_F = @(x) Zn66_Ga66_F(x-3);
plot(Eval,Zn64_Ga64_F(Eval),'g-')
Zn68_Ga68_F = @(x) 1.3 * Zn66_Ga66_F(x+1.5);
plot(Zn68_Ga68_E,Zn68_Ga68_CS,'kd')
plot(Eval,Zn68_Ga68_F(Eval),'k-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
legend('Zn66 data','Zn66 fit','Zn64 data','Zn64 fit','Zn68 data','Zn68 fit')
title('Zn cross sections')


%% Fit secciones eficaces Zn67-Ga67
figure
plot(Zn67_Ga67_E,Zn67_Ga67_CS,'bo'); hold on
Zn67_Ga67_F = fit(Zn67_Ga67_E,Zn67_Ga67_CS,'smoothingspline','SmoothingParam',0.02)
Zn67_Ga67_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,Zn67_Ga67_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Zn67->Ga67 cross sections')
axis([0 300 0 1])

%% Fit secciones eficaces Zn66-Ga65
figure
plot(Zn66_Ga65_E,Zn66_Ga65_CS,'bo'); hold on
Zn66_Ga65_F = fit(Zn66_Ga65_E,Zn66_Ga65_CS,'smoothingspline','SmoothingParam',0.999)
Zn66_Ga65_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,Zn66_Ga65_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Zn66->Ga65 cross sections')
axis([0 300 0 1])

%% Fit secciones eficaces Zn67-Ga66
figure
plot(Zn67_Ga66_E,Zn67_Ga66_CS,'bo'); hold on
Zn67_Ga66_F = fit(Zn67_Ga66_E,Zn67_Ga66_CS,'smoothingspline','SmoothingParam',0.999)
Zn67_Ga66_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,Zn67_Ga66_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Zn67->Ga66 cross sections')
axis([0 300 0 1])

%% Fit secciones eficaces Zn68-Ga67
figure
plot(Zn68_Ga67_E,Zn68_Ga67_CS,'bo'); hold on
Zn68_Ga67_F = fit(Zn68_Ga67_E,Zn68_Ga67_CS,'smoothingspline','SmoothingParam',0.999)
Zn68_Ga67_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,Zn68_Ga67_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Zn68->Ga67 cross sections')
axis([0 300 0 1])

%% Fit secciones eficaces C12_C10
figure
plot(C12_C10_E,C12_C10_CS,'bo'); hold on
C12_C10_F = fit(C12_C10_E,C12_C10_CS,'smoothingspline','SmoothingParam',0.999)
C12_C10_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,C12_C10_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('C12->C10 cross sections')
axis([0 300 0 1])
%axis([0 100 0 0.2]);

%% Fit secciones eficaces  C12_C11
figure
hold off
plot(C12_C11_E,C12_C11_CS,'bo'); hold on
C12_C11_F = fit(C12_C11_E,C12_C11_CS,'smoothingspline','SmoothingParam',0.1)
C12_C11_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,C12_C11_F(Eval),'r-')
axis([0 300 0 0.4])
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('C12->C11 cross sections')

%% Fit secciones eficaces  N14_C11
figure
hold off
plot(N14_C11_E,N14_C11_CS,'bo'); hold on
N14_C11_F = fit(N14_C11_E,N14_C11_CS,'smoothingspline','SmoothingParam',0.999)
N14_C11_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,N14_C11_F(Eval),'r-')
axis([0 200 0 0.4]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('N14->C11 cross sections')

%% Fit secciones eficaces  N14_N13
figure
hold off
plot(N14_N13_E,N14_N13_CS,'bo'); hold on
N14_N13_F = fit(N14_N13_E,N14_N13_CS,'smoothingspline','SmoothingParam',0.2)
N14_N13_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,N14_N13_F(Eval),'r-')
axis([0 200 0 0.1]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('N14->N13 cross sections')

%% Fit secciones eficaces  N14_O14
figure
hold off
plot(N14_O14_E,N14_O14_CS,'bo'); hold on
N14_O14_F = fit(N14_O14_E,N14_O14_CS,'smoothingspline','SmoothingParam',0.99)
N14_O14_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,N14_O14_F(Eval),'r-')
axis([0 200 0 0.1]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('N14->O14 cross sections')

%% Fit secciones eficaces O16_C11
figure
hold off
plot(O16_C11_E,O16_C11_CS,'bo'); hold on
O16_C11_F = fit(O16_C11_E,O16_C11_CS,'smoothingspline','SmoothingParam',0.1)
O16_C11_F.p.coefs(1,:) = [0 0 0 0];
O16_C11_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O16_C11_F(Eval),'r-')
axis([0 200 0 0.1]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('O16->C11 cross sections')

%% Fit secciones eficaces O16_N13
figure
hold off
plot(O16_N13_E,O16_N13_CS,'bo'); hold on
O16_N13_F = fit(O16_N13_E,O16_N13_CS,'smoothingspline','SmoothingParam',0.2)
O16_N13_F.p.coefs(1,:) = [0 0 0 0];
O16_N13_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O16_N13_F(Eval),'r-')
axis([0 100 0 0.2]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('O16->N13 cross sections')


%% Fit secciones eficaces  O16_O15
figure
hold off
plot(O16_O15_E,O16_O15_CS,'bo'); hold on
O16_O15_F = fit(O16_O15_E,O16_O15_CS,'smoothingspline','SmoothingParam',0.2)
O16_O15_F.p.coefs(1,:) = [0 0 0 0];
O16_O15_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O16_O15_F(Eval),'r-')
axis([0 300 0 0.2]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('O16->O15 cross sections')

%% Fit secciones eficaces  O18_F18
figure
hold off
plot(O18_F18_E,O18_F18_CS,'bo'); hold on
O18_F18_F = fit(O18_F18_E,O18_F18_CS,'smoothingspline','SmoothingParam',0.2)
O18_F18_F.p.coefs(1,:) = [0 0 0 0];
O18_F18_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O18_F18_F(Eval),'r-')
axis([0 300 0 1]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('O18->F18 cross sections')

%% Fit secciones eficaces  C44_Sc44
figure
hold off
plot(Ca44_Sc44_E,Ca44_Sc44_CS,'bo'); hold on
Ca44_Sc44_F = fit(Ca44_Sc44_E,Ca44_Sc44_CS,'smoothingspline','SmoothingParam',0.002)
Ca44_Sc44_F.p.coefs(1,:) = [0 0 0 0];
Ca44_Sc44_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Ca44_Sc44_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Ca44_Sc44 cross sections')


%% PG 12C-12C 4.4 MeV
figure
hold off
plot(PG_C12_C12_4_E,PG_C12_C12_4_CS,'bo'); hold on
PG_C12_C12_4_F = fit(PG_C12_C12_4_E,PG_C12_C12_4_CS,'smoothingspline','SmoothingParam',0.2)
PG_C12_C12_4_F.p.coefs(1,:) = [0 0 0 0];
PG_C12_C12_4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,PG_C12_C12_4_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('12C-12C 4.4 MeV cross sections')

%% PG 16O-12C 4.4 MeV
figure
hold off
plot(PG_O16_C12_4_E,PG_O16_C12_4_CS,'bo'); hold on
PG_O16_C12_4_F = fit(PG_O16_C12_4_E,PG_O16_C12_4_CS,'smoothingspline','SmoothingParam',0.2)
PG_O16_C12_4_F.p.coefs(1,:) = [0 0 0 0];
PG_O16_C12_4_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,PG_O16_C12_4_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('16O-12C 4.4 MeV cross sections')

%% PG 14N-14N 1.63 MeV
figure
hold off
plot(PG_N14_N14_1_E,PG_N14_N14_1_CS,'bo'); hold on
PG_N14_N14_1_F = fit(PG_N14_N14_1_E,PG_N14_N14_1_CS,'smoothingspline','SmoothingParam',0.2)
PG_N14_N14_1_F.p.coefs(1,:) = [0 0 0 0];
PG_N14_N14_1_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,PG_N14_N14_1_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('14N-14N 1.63 MeV cross sections')


%% PG 14N-14N 2.31 MeV
figure
hold off
plot(PG_N14_N14_2_E,PG_N14_N14_2_CS,'bo'); hold on
PG_N14_N14_2_F = fit(PG_N14_N14_2_E,PG_N14_N14_2_CS,'smoothingspline','SmoothingParam',0.2)
PG_N14_N14_2_F.p.coefs(1,:) = [0 0 0 0];
PG_N14_N14_2_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,PG_N14_N14_2_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('14N-14N 2.31 MeV cross sections')


%% PG 16 O-16 O 6.13 MeV
figure
hold off
plot(PG_O16_O16_6_E,PG_O16_O16_6_CS,'bo'); hold on
PG_O16_O16_6_F = fit(PG_O16_O16_6_E,PG_O16_O16_6_CS,'smoothingspline','SmoothingParam',0.4)
PG_O16_O16_6_F.p.coefs(1,:) = [0 0 0 0];
PG_O16_O16_6_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,PG_O16_O16_6_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('16 O-16 O 6.13 MeV cross sections')

%% Fit 127I-127Xe MeV
figure
hold off
plot(I127_Xe127_E,I127_Xe127_CS,'bo'); hold on
I127_Xe127_F = fit(I127_Xe127_E,I127_Xe127_CS,'smoothingspline','SmoothingParam',0.9)
I127_Xe127_F.p.coefs(1,:) = [0 0 0 0];
I127_Xe127_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,I127_Xe127_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('I127-Xe127 cross sections')

%% Fit 127I-125Xe MeV
figure
hold off
plot(I127_Xe125_E,I127_Xe125_CS,'bo'); hold on
I127_Xe125_F = fit(I127_Xe125_E,I127_Xe125_CS,'smoothingspline','SmoothingParam',0.9)
I127_Xe125_F.p.coefs(1,:) = [0 0 0 0];
I127_Xe125_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,I127_Xe125_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('I127-Xe125 cross sections')


%% Fit 127I-123Xe MeV
figure
hold off
plot(I127_Xe123_E,I127_Xe123_CS,'bo'); hold on
I127_Xe123_F = fit(I127_Xe123_E,I127_Xe123_CS,'smoothingspline','SmoothingParam',0.9)
I127_Xe123_F.p.coefs(1,:) = [0 0 0 0];
I127_Xe123_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,I127_Xe123_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('I127-Xe123 cross sections')


%% Fit 127I-122Xe MeV
figure
hold off
plot(I127_Xe122_E,I127_Xe122_CS,'bo'); hold on
I127_Xe122_F = fit(I127_Xe122_E,I127_Xe122_CS,'smoothingspline','SmoothingParam',0.99)
I127_Xe122_F.p.coefs(1,:) = [0 0 0 0];
I127_Xe122_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,I127_Xe122_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('I127-Xe122 cross sections')

%% Fit 127I-127mXe MeV
figure
hold off
plot(I127_Xem_E,I127_Xem_CS,'bo'); hold on
I127_Xem_F = fit(I127_Xem_E,I127_Xem_CS,'smoothingspline','SmoothingParam',0.99)
I127_Xem_F.p.coefs(1,:) = [0 0 0 0];
I127_Xem_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,I127_Xem_F(Eval),'r-')
axis([0 300 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('I127-Xe127m cross sections')

%% Fit 23Na-23Mg
figure
hold off
plot(Na23_Mg23_E,Na23_Mg23_CS,'bo'); hold on
Na23_Mg23_F = fit(Na23_Mg23_E,Na23_Mg23_CS,'smoothingspline','SmoothingParam',0.99)
Na23_Mg23_F.p.coefs(1,:) = [0 0 0 0];
Na23_Mg23_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Na23_Mg23_F(Eval),'r-')
axis([0 50 0 1.0]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('Na23Mg23 cross sections')

%% Plot all (water)
figure
hold off
plot(Eval,O16_O15_F(Eval))
hold on
plot(Eval,O16_N13_F(Eval))
plot(Eval,O16_C11_F(Eval))
plot(Eval,Zn64_Ga64_F(Eval))  ; hold on
plot(Eval,Zn66_Ga66_F(Eval))
plot(Eval,Zn68_Ga68_F(Eval))
xlabel('Proton Energy (MeV)')
ylabel('Cross section (barn)');
legend('O16->O15','O16->N13','O16->C11','Zn64->Ga64','Zn66->Ga66','Zn68->Ga68');
title('All cross sections for water + Zn')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis([1 300 1e-4 1]);

%% Plot all (tissue)
figure
hold off
plot(Eval,O16_O15_F(Eval))
hold on
plot(Eval,O16_N13_F(Eval))
plot(Eval,O16_C11_F(Eval))
plot(Eval,N14_N13_F(Eval))
plot(Eval,N14_C11_F(Eval))
plot(Eval,C12_C11_F(Eval))
plot(Eval,C12_C10_F(Eval))
plot(Eval,Zn64_Ga64_F(Eval))
plot(Eval,Zn66_Ga66_F(Eval))
plot(Eval,Zn68_Ga68_F(Eval))
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)');
legend('O16->O15','O16->N13','O16->C11','N14->N13','N14->C11','C12->C11','C12->C10','Zn64->Ga64','Zn66->Ga66','Zn68->Ga68');
title('All cross sections for tissue + Zn')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis([1 300 1e-4 1]);

%% Create fit for stopping power
figure
S_Zn_F = fit(E_keV,S_Zn66,'smoothingspline','SmoothingParam',0.002)
S_Zn_F.p.coefs(1,:) = [0 0 0 0];
S_Zn_F.p.coefs(end,:) = [0 0 0 0];
S_w_F = fit(E_keV,S_w,'smoothingspline','SmoothingParam',0.002)
S_w_F.p.coefs(1,:) = [0 0 0 0];
S_w_F.p.coefs(end,:) = [0 0 0 0];
S_t_F = fit(E_keV,S_tissue,'smoothingspline','SmoothingParam',0.002)
S_t_F.p.coefs(1,:) = [0 0 0 0];
S_t_F.p.coefs(end,:) = [0 0 0 0];
S_b_F = fit(E_keV_bone,S_bone,'smoothingspline','SmoothingParam',0.002)
S_b_F.p.coefs(1,:) = [0 0 0 0];
S_b_F.p.coefs(end,:) = [0 0 0 0];
S_a_F = fit(E_keV_adipose,S_adipose,'smoothingspline','SmoothingParam',0.002)
S_a_F.p.coefs(1,:) = [0 0 0 0];
S_a_F.p.coefs(end,:) = [0 0 0 0];
S_p_F = fit(E_keV_PMMA,S_PMMA,'smoothingspline','SmoothingParam',0.002)
S_p_F.p.coefs(1,:) = [0 0 0 0];
S_p_F.p.coefs(end,:) = [0 0 0 0];
S_pb_F = fit(E_keV_PMMA,S_Pb,'smoothingspline','SmoothingParam',0.002)
S_pb_F.p.coefs(1,:) = [0 0 0 0];
S_pb_F.p.coefs(end,:) = [0 0 0 0];
S_i_F = fit(E_keV_Iodo,S_Iodo,'smoothingspline','SmoothingParam',0.002)
S_i_F.p.coefs(1,:) = [0 0 0 0];
S_i_F.p.coefs(end,:) = [0 0 0 0];
S_CsI_F = fit(E_keV_CsI,S_CsI,'smoothingspline','SmoothingParam',0.002)
S_CsI_F.p.coefs(1,:) = [0 0 0 0];
S_CsI_F.p.coefs(end,:) = [0 0 0 0];
S_NaI_F = fit(E_keV_NaI,S_NaI,'smoothingspline','SmoothingParam',0.002)
S_NaI_F.p.coefs(1,:) = [0 0 0 0];
S_NaI_F.p.coefs(end,:) = [0 0 0 0];
S_Carbon_F = fit(E_keV_carbon,S_Carbon,'smoothingspline','SmoothingParam',0.002)
S_Carbon_F.p.coefs(1,:) = [0 0 0 0];
S_Carbon_F.p.coefs(end,:) = [0 0 0 0];
S_w18_F = fit(E_keV_w18,S_w18,'smoothingspline','SmoothingParam',0.002)
S_w18_F.p.coefs(1,:) = [0 0 0 0];
S_w18_F.p.coefs(end,:) = [0 0 0 0];
loglog(E_keV,S_Zn_F(E_keV),'r-')
hold on;
% loglog(E_keV,S_Zn66,'ro')
% loglog(E_keV,S_w_F(E_keV),'b-')
% loglog(E_keV,S_w,'bo')
% loglog(E_keV,S_t_F(E_keV),'g-')
% loglog(E_keV,S_tissue,'go')
% loglog(E_keV_bone,S_bone,'yo')
% loglog(E_keV_bone,S_b_F(E_keV_bone),'y-')
 loglog(E_keV_adipose,S_adipose,'mo')
 loglog(E_keV_adipose,S_a_F(E_keV_adipose),'m-')
% loglog(E_keV_PMMA,S_PMMA,'co')
% loglog(E_keV_PMMA,S_p_F(E_keV_PMMA),'c-')
%loglog(E_keV_PMMA,S_Pb)
%loglog(E_keV_PMMA,S_pb_F(E_keV_PMMA))
loglog(E_keV_Iodo,S_Iodo,'ko')
loglog(E_keV_Iodo,S_i_F(E_keV_Iodo),'k-')
loglog(E_keV_CsI,S_CsI,'mo')
loglog(E_keV_CsI,S_CsI_F(E_keV_CsI),'m-')
loglog(E_keV_NaI,S_NaI,'bo')
loglog(E_keV_NaI,S_NaI_F(E_keV_NaI),'b-')
loglog(E_keV_w18,S_w18,'yo')
loglog(E_keV_w18,S_w18_F(E_keV_w18),'y-')
legend('Zn66 fit','Zn66 SRIM data', 'Water fit', 'Water SRIM data', 'Tissue fit', 'Tissue SRIM data','Bone fit', 'Bone data', 'Adipose data','Adipose fit','PMMA','fit PMMA', 'Iodo','fit Iodo','CsI','fit CsI','NaI','fir NaI');
xlabel('Proton energy (MeV)');
ylabel('Stopping power (MeV/(cm2/mg))');
%%
close all

