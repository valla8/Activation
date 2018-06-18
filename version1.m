% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es

%% Par?metros a modificar:
clear all
Zn_fraction = 0.1; % Fraccion de Zn (tanto por uno)
dx = 0.01; % Espaciado de la malla (cm)
E0 = 58; % Energ?a inicial del haz (MeV)

%% Cargar datos
% Secciones de EXFOR.
load('CrossSections.mat');
%C12_C10_CS=C12_C10_CS./1000;
C12_C11_CS=C12_C11_CS./1000;
PG_C12_C12_4_CS=PG_C12_C12_4_CS'./1000;
PG_O16_C12_4_CS=PG_O16_C12_4_CS'./1000;
PG_N14_N14_1_CS=PG_N14_N14_1_CS'./1000;
PG_N14_N14_2_CS=PG_N14_N14_2_CS'./1000;
PG_O16_O16_6_CS=PG_O16_O16_6_CS'./1000;
PG_C12_C12_4_E=PG_C12_C12_4_E';
PG_O16_C12_4_E=PG_O16_C12_4_E';
PG_N14_N14_1_E=PG_N14_N14_1_E';
PG_N14_N14_2_E=PG_N14_N14_2_E';
PG_O16_O16_6_E=PG_O16_O16_6_E';
% Vidas medias en s
load('MeanLives.mat');
landa_C10 = log(2) / T_C10;
landa_C11 = log(2) / T_C11;
landa_Ga64 = log(2) / T_Ga64;
landa_Ga66 = log(2) / T_Ga66;
landa_Ga68 = log(2) / T_Ga68;
landa_N13 = log(2) / T_N13;
landa_O15 = log(2) / T_O15;
landa_Sc44 = log(2) / T_Sc44;

% Natural abundances
O16_ab = 0.99729; %1
N14_ab = 0.996; %2
C12_ab = 0.989; %3
Zn64_ab = 0.492; %4
Zn66_ab = 0.277; %5
Zn68_ab = 0.185; %6
Ca44_ab = 0.0223233841; %
ab = [O16_ab N14_ab C12_ab Zn64_ab Zn66_ab Zn68_ab];

% Tissue compositions (fuente: Zhu 2011, ICRU Report #63, ICRU Report #44)
% Material H (%) C (%) N (%) O (%) Zn64(%) Zn66(%) Zn68(%)  Density (g ml?1)
% Soft tissue 0.102 0.143 0.034 0.708 0 0 0 1.06
% Water: 0.0305 0 0 0.9695 0 0 0 1.00
Comp_tissue = [0.619977 0.118053 0.0205881 0.240345 0 0 ];
Comp_tissue_2 = [0.619977 0.118053 0.0205881 0.240345 0 0 0];
Comp_water = [0.0305 0 0 0.9695 0 0 0];
Comp_Zn = [0 0 0 0 Zn64_ab Zn66_ab Zn68_ab];
Comp_tissue_Zn =Zn_fraction*Comp_Zn + (1-Zn_fraction)* Comp_tissue_2;!
Comp_water_Zn = Zn_fraction*Comp_Zn + (1-Zn_fraction)*Comp_water;
% Material H (%) C (%) N (%) O (%) P (%) Ca(%)  Density (g ml?1)
Comp_bone=[0.475402    0.121899    0.0304112    0.282845    0.0343791    0.0531337];
W_ele = [1.00794 12.011 14.00674 15.9994 30.973762 40.078];
Comp_adipose = [0.634657 0.284  0.0030464 0.0777462 0 0];
Comp_PMMA= [0.535002   0.332244 0  0.132754 0 0];

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

%% Fit secciones eficaces C12_C10
figure
plot(C12_C10_E,C12_C10_CS,'bo'); hold on
C12_C10_F = fit(C12_C10_E,C12_C10_CS,'smoothingspline','SmoothingParam',0.999)
C12_C10_F.p.coefs(1,:) = [0 0 0 0]
plot(Eval,C12_C10_F(Eval),'r-')
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('C12->C10 cross sections')
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
N14_N13_F = fit(N14_N13_E,N14_N13_CS,'smoothingspline','SmoothingParam',0.5)
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
O16_N13_F = fit(O16_N13_E,O16_N13_CS,'smoothingspline','SmoothingParam',0.999)
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
O16_O15_F = fit(O16_O15_E,O16_O15_CS,'smoothingspline','SmoothingParam',0.002)
O16_O15_F.p.coefs(1,:) = [0 0 0 0];
O16_O15_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,O16_O15_F(Eval),'r-')
axis([0 300 0 0.2]);
xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
title('O16->O15 cross sections')

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
loglog(E_keV,S_Zn_F(E_keV),'r-')
hold on;
loglog(E_keV,S_Zn66,'ro')
loglog(E_keV,S_w_F(E_keV),'b-')
loglog(E_keV,S_w,'bo')
loglog(E_keV,S_t_F(E_keV),'g-')
loglog(E_keV,S_tissue,'go')
loglog(E_keV_bone,S_bone,'yo')
loglog(E_keV_bone,S_b_F(E_keV_bone),'y-')
loglog(E_keV_adipose,S_adipose,'mo')
loglog(E_keV_adipose,S_a_F(E_keV_adipose),'m-')
loglog(E_keV_PMMA,S_PMMA,'co')
loglog(E_keV_PMMA,S_p_F(E_keV_PMMA),'c-')
loglog(E_keV_PMMA,S_Pb)
loglog(E_keV_PMMA,S_pb_F(E_keV_PMMA))
legend('Zn66 fit','Zn66 SRIM data', 'Water fit', 'Water SRIM data', 'Tissue fit', 'Tissue SRIM data','Bone fit', 'Bone data', 'Adipose data','Adipose fit','PMMA','fit PMMA', 'Pb');
xlabel('Proton energy (MeV)');
ylabel('Stopping power (MeV/(cm2/mg))');
%%
close all

%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
PMMA_Weight=100.12 %g/mol
rho_w = 1; % g/cm3
rho_Zn = 7.14; %g/cm3
rho_tissue = 1.1; %g/cm3
rho_bone = 1.85; %g/cm3
rho_adipose = 0.92; %g/cm3
rho_PMMA= 1.18; %g/cm3
!rho_tissue_A = 4.6243E+22; % atoms/cm3
rho_tissue_A = AvNmbr*rho_tissue/sum(Comp_tissue.*W_ele);  % atoms/cm3
rho_bone_A = AvNmbr*rho_bone/sum(Comp_bone.*W_ele);  % atoms/cm3
rho_adipose_A = AvNmbr*rho_adipose/sum(Comp_adipose.*W_ele);  % atoms/cm3
rho_PMMA_A = AvNmbr*rho_PMMA/PMMA_Weight;  % atoms/cm3
rho_w_A = (1-Zn_fraction) * rho_w * AvNmbr / waterMolecularWeight; % molecules / cm3
ZnAtomicWeight = 65.38; % g/mol
rho_Zn_A = Zn_fraction * rho_Zn * AvNmbr / ZnAtomicWeight; % molecules / cm3
rho_O16_A = rho_w_A * O16_ab; % atoms/cm3
rho_O16_Ab = rho_bone_A * Comp_bone(4) * O16_ab;
rho_N14_Ab = rho_bone_A * Comp_bone(3) * N14_ab;
rho_C12_Ab = rho_bone_A * Comp_bone(2) * C12_ab;
rho_Ca44_Ab = rho_bone_A * Comp_bone(6) * Ca44_ab;
rho_O16_At = rho_tissue_A * Comp_tissue_Zn(4) * O16_ab;
rho_N14_At = rho_tissue_A * Comp_tissue_Zn(3) * N14_ab;
rho_C12_At = rho_tissue_A * Comp_tissue_Zn(2) * C12_ab;
rho_O16_Ap = rho_PMMA_A * Comp_PMMA(4) * O16_ab;
rho_N14_Ap = rho_PMMA_A * Comp_PMMA(3) * N14_ab;
rho_C12_Ap = rho_PMMA_A * Comp_PMMA(2) * C12_ab;
rho_O16_Aa = rho_adipose_A * Comp_adipose(4) * O16_ab;
rho_N14_Aa = rho_adipose_A * Comp_adipose(3) * N14_ab;
rho_C12_Aa = rho_adipose_A * Comp_adipose(2) * C12_ab;
pps=1; %proton per second

x = 0:dx:5; % posiciones en cm.
E = nan(size(x));
Et = nan(size(x));
Eb = nan(size(x));
Ea = nan(size(x));
Ep = nan(size(x));

Ddep = zeros(size(E));
Ddept = zeros(size(E));
Ddepb = zeros(size(E));
Ddepa = zeros(size(E));
Ddepp = zeros(size(E));

currentE = E0;
currentEt = E0;
currentEb = E0;
currentEa = E0;
currentEp = E0;


% In water (full + simplified versions)
Y64 = nan(size(x));
Y66 = nan(size(x));
Y68 = nan(size(x));
Y64s = nan(size(x));
Y66s = nan(size(x));
Y68s = nan(size(x));
Y_O16_C11 = nan(size(x));
Y_O16_N13 = nan(size(x));
Y_O16_O15 = nan(size(x));
Y_O16_C11s = nan(size(x));
Y_O16_N13s = nan(size(x));
Y_O16_O15s = nan(size(x));

% In tissue (simplified form only)
Y64t = nan(size(x));
Y66t = nan(size(x));
Y68t = nan(size(x));
Y_O16_C11t = zeros(size(x));
Y_O16_N13t = zeros(size(x));
Y_O16_O15t = zeros(size(x));
Y_N14_N13t = zeros(size(x));
Y_N14_O14t = zeros(size(x));
Y_N14_C11t = zeros(size(x));   
Y_C12_C11t = zeros(size(x));     
Y_C12_C10t = zeros(size(x)); 
Y_PG_C12_C12_4t = zeros(size(x));
Y_PG_O16_C12_4t = zeros(size(x));
Y_PG_N14_N14_1t = zeros(size(x));
Y_PG_N14_N14_2t = zeros(size(x));
Y_PG_O16_O16_6t = zeros(size(x));

% In bone (simplified form only)
Y_O16_C11b = zeros(size(x));
Y_O16_N13b = zeros(size(x));
Y_O16_O15b = zeros(size(x));
Y_N14_N13b = zeros(size(x));
Y_N14_O14b = zeros(size(x));
Y_N14_C11b = zeros(size(x));   
Y_C12_C11b = zeros(size(x));     
Y_C12_C10b = zeros(size(x)); 
Y_Ca44_Sc44 = zeros(size(x)); 
Y_PG_C12_C12_4b = zeros(size(x));
Y_PG_O16_C12_4b = zeros(size(x));
Y_PG_N14_N14_1b = zeros(size(x));
Y_PG_N14_N14_2b = zeros(size(x));
Y_PG_O16_O16_6b = zeros(size(x));

% In adipose (simplified form only)
Y_O16_C11a = zeros(size(x));
Y_O16_N13a = zeros(size(x));
Y_O16_O15a = zeros(size(x));
Y_N14_N13a = zeros(size(x));
Y_N14_O14a = zeros(size(x));
Y_N14_C11a = zeros(size(x));   
Y_C12_C11a = zeros(size(x));     
Y_C12_C10a = zeros(size(x)); 
Y_PG_C12_C12_4a = zeros(size(x));
Y_PG_O16_C12_4a = zeros(size(x));
Y_PG_N14_N14_1a = zeros(size(x));
Y_PG_N14_N14_2a = zeros(size(x));
Y_PG_O16_O16_6a = zeros(size(x));

% In PMMA (simplified form only)
Y_O16_C11p = zeros(size(x));
Y_O16_N13p = zeros(size(x));
Y_O16_O15p = zeros(size(x));
Y_N14_N13p = zeros(size(x));
Y_N14_O14p = zeros(size(x));
Y_N14_C11p = zeros(size(x));   
Y_C12_C11p = zeros(size(x));     
Y_C12_C10p = zeros(size(x));
Y_PG_C12_C12_4p = zeros(size(x));
Y_PG_O16_C12_4p = zeros(size(x));
Y_PG_N14_N14_1p = zeros(size(x));
Y_PG_N14_N14_2p = zeros(size(x));
Y_PG_O16_O16_6p = zeros(size(x));

for i=1:(numel(x)-1)
    S_w = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Zn = max(0,1000*S_Zn_F(currentE*1000)); % MeV/(g/cm2) 
    S_Znt = max(0,1000*S_Zn_F(currentEt*1000)); % MeV/(g/cm2) 
    S_t = max(0,1000*S_t_F(currentEt*1000)); % MeV/(g/cm2)
    S_b = max(0,1000*S_b_F(currentEb*1000));
    S_a = max(0,1000*S_a_F(currentEa*1000));
    S_p = max(0,1000*S_p_F(currentEp*1000));
    
    S1 = (S_w*rho_w*(1-Zn_fraction) + S_Zn*rho_Zn*Zn_fraction); % MeV/cm
    S1t = (S_t*rho_tissue); %*(1-Zn_fraction) + S_Znt*rho_Zn*Zn_fraction); % MeV/cm
    S1b = (S_b*rho_bone); % MeV/cm
    S1a = (S_a*rho_adipose); % MeV/cm
    S1p = (S_p*rho_PMMA); % MeV/cm
    
    deltaE = dx*S1; % MeV
    deltaEt = dx*S1t; % MeV
    deltaEb = dx*S1b; % MeV
    deltaEa = dx*S1a; %MeV
    deltaEp = dx*S1p; %MeV
    
    E(i) = currentE; % MeV
    Et(i) = currentEt; % MeV
    Eb(i) = currentEb;  % MeV
    Ea(i) = currentEa; % MeV
    Ep(i) = currentEp; % MeV
    currentE = currentE - deltaE; % MeV
    currentEt = currentEt - deltaEt; % MeV
    currentEb = currentEb - deltaEb; % MeV
    currentEa = currentEa - deltaEa; %MeV
    currentEp = currentEp - deltaEp; %MeV
    
    S_w2 = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Zn2 = max(0,1000*S_Zn_F(currentE*1000)); % MeV/(g/cm2)
    S_t2 = max(0,1000*S_t_F(currentEt*1000)); % MeV/(g/cm2)
    S_Zn2t = max(0,1000*S_Zn_F(currentEt*1000)); % MeV/(g/cm2)
    S_b2 = max(0,1000*S_b_F(currentEb*1000));% MeV/(g/cm2)
    S_a2 = max(0,1000*S_a_F(currentEa*1000));% MeV/(g/cm2)
    S_p2 = max(0,1000*S_p_F(currentEp*1000));% MeV/(g/cm2)
    
    S2 = (S_w2*rho_w*(1-Zn_fraction) + S_Zn2*rho_Zn*Zn_fraction); % MeV/cm
    S2t = (S_t2*rho_tissue); %*(1-Zn_fraction) + S_Zn2t*rho_Zn*Zn_fraction); % MeV/cm
    S2b = (S_b2*rho_bone); % MeV/cm3
    S2a = (S_a2*rho_adipose); % MeV/cm3
    S2p = (S_p2*rho_PMMA); % MeV/cm3
    
    Stt=(S2t+S1t)/2;
    Stb=(S2b+S1b)/2;
    Sta=(S2a+S1a)/2;
    Stp=(S2p+S1p)/2;
    
    Ddep(i) = deltaE; % MeV
    Ddept(i) = dx*Stt; % MeV
    Ddepb(i) = dx*Stb;  % MeV
    Ddepa(i) = dx*Sta; %MeV
    Ddepp(i) = dx*Stp; %MeV
    
    % Zn part
    E1 = E(i);
    E1t = Et(i);
    E1b = Eb(i);
    E1a = Ea(i);
    E1p = Ep(i);
    E2 = currentE;
    E2t = currentEt;
    E2a = currentEa;
    E2b = currentEb;
    E2p = currentEp;
    E2E1 = [currentE E(i)];
    E2E1t = [currentEt Et(i)];
    E2E1b = [currentEb Eb(i)];
    
    % Water (full + simplified)
    Y64(i) = rho_Zn_A * Zn64_ab * trapz(E2E1, [Zn64_Ga64_F(E1)*1e-24/S1 Zn64_Ga64_F(E2)*1e-24/S2]);
    Y66(i) = rho_Zn_A * Zn66_ab * trapz(E2E1, [Zn66_Ga66_F(E1)*1e-24/S1 Zn66_Ga66_F(E2)*1e-24/S2]);
    Y68(i) = rho_Zn_A * Zn68_ab * trapz(E2E1, [Zn68_Ga68_F(E1)*1e-24/S1 Zn68_Ga68_F(E2)*1e-24/S2]);
    Y_O16_C11(i) = rho_O16_A * trapz(E2E1, [O16_C11_F(E1)*1e-24/S1 O16_C11_F(E2)*1e-24/S2]);
    Y_O16_N13(i) = rho_O16_A * trapz(E2E1, [O16_N13_F(E1)*1e-24/S1 O16_N13_F(E2)*1e-24/S2]);
    Y_O16_O15(i) = rho_O16_A * trapz(E2E1, [O16_O15_F(E1)*1e-24/S1 O16_O15_F(E2)*1e-24/S2]);
    sigma_64_mean = 0.5 * (Zn64_Ga64_F(E1) + Zn64_Ga64_F(E2));
    sigma_66_mean = 0.5 * (Zn66_Ga66_F(E1) + Zn66_Ga66_F(E2));
    sigma_68_mean = 0.5 * (Zn68_Ga68_F(E1) + Zn68_Ga68_F(E2));
    sigma_C11_mean = 0.5 * (O16_C11_F(E1) + O16_C11_F(E2));
    sigma_N13_mean = 0.5 * (O16_N13_F(E1) + O16_N13_F(E2));
    sigma_O15_mean = 0.5 * (O16_O15_F(E1) + O16_O15_F(E2));
    Y64s(i) = rho_Zn_A * Zn64_ab * sigma_64_mean * 1e-24 * dx;
    Y66s(i) = rho_Zn_A * Zn66_ab * sigma_66_mean * 1e-24 * dx;
    Y68s(i) = rho_Zn_A * Zn68_ab * sigma_68_mean * 1e-24 * dx;
    Y_O16_C11s(i) = rho_O16_A * sigma_C11_mean * 1e-24 * dx;
    Y_O16_N13s(i) = rho_O16_A * sigma_N13_mean * 1e-24 * dx;
    Y_O16_O15s(i) = rho_O16_A * sigma_O15_mean * 1e-24 * dx;
    
    % Tissue (simplified only)
    sigma_64_meant = 0.5 * (max(0,Zn64_Ga64_F(E1t)) + max(0,Zn64_Ga64_F(E2t)));
    sigma_66_meant = 0.5 * (max(0,Zn66_Ga66_F(E1t)) + max(0,Zn66_Ga66_F(E2t)));
    sigma_68_meant = 0.5 * (max(0,Zn68_Ga68_F(E1t)) + max(0,Zn68_Ga68_F(E2t)));
    sigma_C11_meant = 0.5 * (max(0,O16_C11_F(E1t)) + max(0,O16_C11_F(E2t)));
    sigma_N13_meant = 0.5 * (max(0,O16_N13_F(E1t)) + max(0,O16_N13_F(E2t)));
    sigma_O15_meant = 0.5 * (max(0,O16_O15_F(E1t)) + max(0,O16_O15_F(E2t)));
    sigma_N14_N13t = 0.5 * (max(0,N14_N13_F(E1t)) + max(0,N14_N13_F(E2t)));
    sigma_N14_O14t = 0.5 * (max(0,N14_O14_F(E1t)) + max(0,N14_O14_F(E2t)));
    sigma_N14_C11t = 0.5 * (max(0,N14_C11_F(E1t)) + max(0,N14_C11_F(E2t)));
    sigma_C12_C11t = 0.5 * (max(0,C12_C11_F(E1t)) + max(0,C12_C11_F(E2t)));
    sigma_C12_C10t = 0.5 * (max(0,C12_C10_F(E1t)) + max(0,C12_C10_F(E2t)));
    sigma_PG_C12_4t = 0.5 * (max(0,PG_C12_C12_4_F(E1t)) + max(0,PG_C12_C12_4_F(E2t)));
    sigma_PG_O16_4t = 0.5 * (max(0,PG_O16_C12_4_F(E1t)) + max(0,PG_O16_C12_4_F(E2t)));
    sigma_PG_N14_1t = 0.5 * (max(0,PG_N14_N14_1_F(E1t)) + max(0,PG_N14_N14_1_F(E2t)));
    sigma_PG_N14_2t = 0.5 * (max(0,PG_N14_N14_2_F(E1t)) + max(0,PG_N14_N14_2_F(E2t)));
    sigma_PG_O16_6t = 0.5 * (max(0,PG_O16_O16_6_F(E1t)) + max(0,PG_O16_O16_6_F(E2t)));
    Y64t(i) = rho_Zn_A * Zn64_ab * sigma_64_meant * 1e-24 * dx;
    Y66t(i) = rho_Zn_A * Zn66_ab * sigma_66_meant * 1e-24 * dx;
    Y68t(i) = rho_Zn_A * Zn68_ab * sigma_68_meant * 1e-24 * dx;
    Y_O16_C11t(i) = pps * rho_O16_At * sigma_C11_meant * 1e-24 * dx;
    Y_O16_N13t(i) = pps * rho_O16_At * sigma_N13_meant * 1e-24 * dx;
    Y_O16_O15t(i) = pps * rho_O16_At * sigma_O15_meant * 1e-24 * dx;
    Y_N14_N13t(i) = pps * rho_N14_At * sigma_N14_N13t * 1e-24 * dx;
    Y_N14_C11t(i) = pps * rho_N14_At * sigma_N14_C11t * 1e-24 * dx;! 
    Y_N14_O14t(i) = pps * rho_N14_At * sigma_N14_O14t * 1e-24 * dx;!
    Y_C12_C11t(i) = pps * rho_C12_At * sigma_C12_C11t * 1e-24 * dx;
    Y_C12_C10t(i) = pps * rho_C12_At * sigma_C12_C10t * 1e-24 * dx;
    Y_PG_C12_C12_4t(i) = pps * rho_C12_At * sigma_PG_C12_4t * 1e-24 * dx;
    Y_PG_O16_C12_4t(i) = pps * rho_O16_At * sigma_PG_O16_4t * 1e-24 * dx;
    Y_PG_N14_N14_1t(i) = pps * rho_N14_At * sigma_PG_N14_1t * 1e-24 * dx;
    Y_PG_N14_N14_2t(i) = pps * rho_N14_At * sigma_PG_N14_2t * 1e-24 * dx;
    Y_PG_O16_O16_6t(i) = pps * rho_O16_At * sigma_PG_O16_6t * 1e-24 * dx;
    
    
    % Bone (simplified only)

    sigma_C11_meanb = 0.5 * (max(0,O16_C11_F(E1b)) + max(0,O16_C11_F(E2b)));
    sigma_N13_meanb = 0.5 * (max(0,O16_N13_F(E1b)) + max(0,O16_N13_F(E2b)));
    sigma_O15_meanb = 0.5 * (max(0,O16_O15_F(E1b)) + max(0,O16_O15_F(E2b)));
    sigma_N14_N13b = 0.5 * (max(0,N14_N13_F(E1b)) + max(0,N14_N13_F(E2b)));
    sigma_N14_O14b = 0.5 * (max(0,N14_O14_F(E1b)) + max(0,N14_O14_F(E2b)));
    sigma_N14_C11b = 0.5 * (max(0,N14_C11_F(E1b)) + max(0,N14_C11_F(E2b)));
    sigma_C12_C11b = 0.5 * (max(0,C12_C11_F(E1b)) + max(0,C12_C11_F(E2b)));
    sigma_C12_C10b = 0.5 * (max(0,C12_C10_F(E1b)) + max(0,C12_C10_F(E2b)));
    sigma_Ca44_Sc44b = 0.5 * (max(0,Ca44_Sc44_F(E1b)) + max(0,Ca44_Sc44_F(E2b)));
    sigma_PG_C12_4b = 0.5 * (max(0,PG_C12_C12_4_F(E1b)) + max(0,PG_C12_C12_4_F(E2b)));
    sigma_PG_O16_4b = 0.5 * (max(0,PG_O16_C12_4_F(E1b)) + max(0,PG_O16_C12_4_F(E2b)));
    sigma_PG_N14_1b = 0.5 * (max(0,PG_N14_N14_1_F(E1b)) + max(0,PG_N14_N14_1_F(E2b)));
    sigma_PG_N14_2b = 0.5 * (max(0,PG_N14_N14_2_F(E1b)) + max(0,PG_N14_N14_2_F(E2b)));
    sigma_PG_O16_6b = 0.5 * (max(0,PG_O16_O16_6_F(E1b)) + max(0,PG_O16_O16_6_F(E2b)));
    Y_O16_C11b(i) = pps * rho_O16_Ab * sigma_C11_meanb * 1e-24 * dx;
    Y_O16_N13b(i) = pps * rho_O16_Ab * sigma_N13_meanb * 1e-24 * dx;
    Y_O16_O15b(i) = pps * rho_O16_Ab * sigma_O15_meanb * 1e-24 * dx;
    Y_N14_N13b(i) = pps * rho_N14_Ab * sigma_N14_N13b * 1e-24 * dx;
    Y_N14_C11b(i) = pps * rho_N14_Ab * sigma_N14_C11b * 1e-24 * dx;! 
    Y_N14_O14b(i) = pps * rho_N14_Ab * sigma_N14_O14b* 1e-24 * dx;!
    Y_C12_C11b(i) = pps * rho_C12_Ab * sigma_C12_C11b * 1e-24 * dx;
    Y_C12_C10b(i) = pps * rho_C12_Ab * sigma_C12_C10b * 1e-24 * dx;
    Y_Ca44_Sc44(i) = pps * rho_Ca44_Ab * sigma_Ca44_Sc44b * 1e-24 * dx;
    Y_PG_C12_C12_4b(i) = pps * rho_C12_Ab * sigma_PG_C12_4b * 1e-24 * dx;
    Y_PG_O16_C12_4b(i) = pps * rho_O16_Ab * sigma_PG_O16_4b * 1e-24 * dx;
    Y_PG_N14_N14_1b(i) = pps * rho_N14_Ab * sigma_PG_N14_1b * 1e-24 * dx;
    Y_PG_N14_N14_2b(i) = pps * rho_N14_Ab * sigma_PG_N14_2b * 1e-24 * dx;
    Y_PG_O16_O16_6b(i) = pps * rho_O16_Ab * sigma_PG_O16_6b * 1e-24 * dx;
    
        % Adipose (simplified only)

    sigma_C11_meana = 0.5 * (max(0,O16_C11_F(E1a)) + max(0,O16_C11_F(E2a)));
    sigma_N13_meana = 0.5 * (max(0,O16_N13_F(E1a)) + max(0,O16_N13_F(E2a)));
    sigma_O15_meana = 0.5 * (max(0,O16_O15_F(E1a)) + max(0,O16_O15_F(E2a)));
    sigma_N14_N13a = 0.5 * (max(0,N14_N13_F(E1a)) + max(0,N14_N13_F(E2a)));
    sigma_N14_O14a = 0.5 * (max(0,N14_O14_F(E1a)) + max(0,N14_O14_F(E2a)));
    sigma_N14_C11a = 0.5 * (max(0,N14_C11_F(E1a)) + max(0,N14_C11_F(E2a)));
    sigma_C12_C11a = 0.5 * (max(0,C12_C11_F(E1a)) + max(0,C12_C11_F(E2a)));
    sigma_C12_C10a = 0.5 * (max(0,C12_C10_F(E1a)) + max(0,C12_C10_F(E2a)));
    sigma_PG_C12_4a = 0.5 * (max(0,PG_C12_C12_4_F(E1a)) + max(0,PG_C12_C12_4_F(E2a)));
    sigma_PG_O16_4a = 0.5 * (max(0,PG_O16_C12_4_F(E1a)) + max(0,PG_O16_C12_4_F(E2a)));
    sigma_PG_N14_1a = 0.5 * (max(0,PG_N14_N14_1_F(E1a)) + max(0,PG_N14_N14_1_F(E2a)));
    sigma_PG_N14_2a = 0.5 * (max(0,PG_N14_N14_2_F(E1a)) + max(0,PG_N14_N14_2_F(E2a)));
    sigma_PG_O16_6a = 0.5 * (max(0,PG_O16_O16_6_F(E1a)) + max(0,PG_O16_O16_6_F(E2a)));
    Y_O16_C11a(i) = pps * rho_O16_Aa * sigma_C11_meana * 1e-24 * dx;
    Y_O16_N13a(i) = pps * rho_O16_Aa * sigma_N13_meana * 1e-24 * dx;
    Y_O16_O15a(i) = pps * rho_O16_Aa * sigma_O15_meana * 1e-24 * dx;
    Y_N14_N13a(i) = pps * rho_N14_Aa * sigma_N14_N13a * 1e-24 * dx;
    Y_N14_C11a(i) = pps * rho_N14_Aa * sigma_N14_C11a * 1e-24 * dx;! 
    Y_N14_O14a(i) = pps * rho_N14_Aa * sigma_N14_O14a* 1e-24 * dx;!
    Y_C12_C11a(i) = pps * rho_C12_Aa * sigma_C12_C11a * 1e-24 * dx;
    Y_C12_C10a(i) = pps * rho_C12_Aa * sigma_C12_C10a * 1e-24 * dx;
    Y_PG_C12_C12_4a(i) = pps * rho_C12_Aa * sigma_PG_C12_4a * 1e-24 * dx;
    Y_PG_O16_C12_4a(i) = pps * rho_O16_Aa * sigma_PG_O16_4a * 1e-24 * dx;
    Y_PG_N14_N14_1a(i) = pps * rho_N14_Aa * sigma_PG_N14_1a * 1e-24 * dx;
    Y_PG_N14_N14_2a(i) = pps * rho_N14_Aa * sigma_PG_N14_2a * 1e-24 * dx;
    Y_PG_O16_O16_6a(i) = pps * rho_O16_Aa * sigma_PG_O16_6a * 1e-24 * dx;
    
    
            % PMMA (simplified only)

    sigma_C11_meanp = 0.5 * (max(0,O16_C11_F(E1p)) + max(0,O16_C11_F(E2p)));
    sigma_N13_meanp = 0.5 * (max(0,O16_N13_F(E1p)) + max(0,O16_N13_F(E2p)));
    sigma_O15_meanp = 0.5 * (max(0,O16_O15_F(E1p)) + max(0,O16_O15_F(E2p)));
    sigma_N14_N13p = 0.5 * (max(0,N14_N13_F(E1p)) + max(0,N14_N13_F(E2p)));
    sigma_N14_O14p = 0.5 * (max(0,N14_O14_F(E1p)) + max(0,N14_O14_F(E2p)));
    sigma_N14_C11p = 0.5 * (max(0,N14_C11_F(E1p)) + max(0,N14_C11_F(E2p)));
    sigma_C12_C11p = 0.5 * (max(0,C12_C11_F(E1p)) + max(0,C12_C11_F(E2p)));
    sigma_C12_C10p = 0.5 * (max(0,C12_C10_F(E1p)) + max(0,C12_C10_F(E2p)));
    sigma_PG_C12_4p = 0.5 * (max(0,PG_C12_C12_4_F(E1p)) + max(0,PG_C12_C12_4_F(E2p)));
    sigma_PG_O16_4p = 0.5 * (max(0,PG_O16_C12_4_F(E1p)) + max(0,PG_O16_C12_4_F(E2p)));
    sigma_PG_N14_1p = 0.5 * (max(0,PG_N14_N14_1_F(E1p)) + max(0,PG_N14_N14_1_F(E2p)));
    sigma_PG_N14_2p = 0.5 * (max(0,PG_N14_N14_2_F(E1p)) + max(0,PG_N14_N14_2_F(E2p)));
    sigma_PG_O16_6p = 0.5 * (max(0,PG_O16_O16_6_F(E1p)) + max(0,PG_O16_O16_6_F(E2p)));
    Y_O16_C11p(i) = pps * rho_O16_Ap * sigma_C11_meanp * 1e-24 * dx;
    Y_O16_N13p(i) = pps * rho_O16_Ap * sigma_N13_meanp * 1e-24 * dx;
    Y_O16_O15p(i) = pps * rho_O16_Ap * sigma_O15_meanp * 1e-24 * dx;
    Y_N14_N13p(i) = pps * rho_N14_Ap * sigma_N14_N13p* 1e-24 * dx;
    Y_N14_C11p(i) = pps * rho_N14_Ap * sigma_N14_C11p * 1e-24 * dx;! 
    Y_N14_O14p(i) = pps * rho_N14_Ap * sigma_N14_O14p* 1e-24 * dx;!
    Y_C12_C11p(i) = pps * rho_C12_Ap * sigma_C12_C11p * 1e-24 * dx;
    Y_C12_C10p(i) = pps * rho_C12_Ap * sigma_C12_C10p * 1e-24 * dx;
    Y_PG_C12_C12_4p(i) = pps * rho_C12_Ap * sigma_PG_C12_4p * 1e-24 * dx;
    Y_PG_O16_C12_4p(i) = pps * rho_O16_Ap * sigma_PG_O16_4p * 1e-24 * dx;
    Y_PG_N14_N14_1p(i) = pps * rho_N14_Ap * sigma_PG_N14_1p * 1e-24 * dx;
    Y_PG_N14_N14_2p(i) = pps * rho_N14_Ap * sigma_PG_N14_2p * 1e-24 * dx;
    Y_PG_O16_O16_6p(i) = pps * rho_O16_Ap * sigma_PG_O16_6p * 1e-24 * dx;
    
    
end

%% Create plots

% Figure in water
figure
%subplot(2,1,1)
%plot(x,E)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Water + Zn');
%hold on
plot(x,100*Ddep)
legend('Dose')
set(gca,'FontSize',14)
axis([0 30 0 (max(100*Ddep)+50)]);
%subplot(2,1,2)
yyaxis left
title('Yields of different species (per incoming proton)');
hold on
plot(x,Y64,'b-'); hold on
plot(x,Y66,'r-')
plot(x,Y68,'g-')
plot(x,Y_O16_C11,'k'); hold on
plot(x,Y_O16_N13,'c')
plot(x,Y_O16_O15,'m')
legend('Zn64','Zn66','Zn68','C11','N13','O15','Location', 'Southeast');
ylabel('Yield');
plot(x,Y64s,'b:'); hold on
plot(x,Y66s,'r:')
plot(x,Y68s,'g:')
plot(x,Y_O16_C11s,'k:'); hold on
plot(x,Y_O16_N13s,'c:')
plot(x,Y_O16_O15s,'m:')
set(gca,'FontSize',14)
[f,g]=min(Ddep);
axis([0  (ceil(g*dx)) 0 0.3e-4]);

% Figure in tisssue
figure
%subplot(2,1,1)
%plot(x,Et)
yyaxis right
xlabel('Depth (cm)');
ylabel(' Dose (a.u.)')
title('Tissue');
hold on
plot(x,100*Ddept)
legend('Dose')
set(gca,'FontSize',14)
axis([0 17 0 (max(100*Ddept)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11t = Y_O16_C11t + Y_C12_C11t + Y_N14_C11t;
Y_N13t = Y_O16_N13t + Y_N14_N13t;
Y_O15t = Y_O16_O15t;
Y_C10t = Y_C12_C10t;
plot(x,Y_C11t,'k'); hold on
plot(x,Y_N13t,'c')
plot(x,Y_O15t,'m')
plot(x,Y_C10t,'y');
legend('C11','N13','O15','C10','Location', 'northwest');
[f,g]=min(Ddept);
axis([0 (ceil(g*dx)) 0 max(max(Y_C11t,Y_O15t)+0.2*max(Y_C11t,Y_O15t))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

% Figure in adipose
figure
%subplot(2,1,1)
%plot(x,Ea)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Adipose');
hold on
plot(x,100*Ddepa)
legend('Dose')
set(gca,'FontSize',14)
axis([0 20 0 (max(100*Ddepa)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11a = Y_O16_C11a + Y_C12_C11a + Y_N14_C11a;
Y_N13a = Y_O16_N13a + Y_N14_N13a;
Y_O15a = Y_O16_O15a;
Y_C10a = Y_C12_C10a;
plot(x,Y_C11a,'k'); hold on
plot(x,Y_N13a,'c')
plot(x,Y_O15a,'m')
plot(x,Y_C10a,'y');
legend('C11','N13','O15','C10','Location', 'northwest');
[f,g]=min(Ddepa);
axis([0 (ceil(g*dx)) 0 max(max(Y_C11a,Y_O15a)+0.2*max(Y_C11a,Y_O15a))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

% Figure in bone
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Bone');
hold on
plot(x,100*Ddepb)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(100*Ddepb)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11b = Y_O16_C11b + Y_C12_C11b + Y_N14_C11b;
Y_N13b = Y_O16_N13b + Y_N14_N13b;
Y_O15b = Y_O16_O15b;
Y_C10b = Y_C12_C10b;
Y_Sc44b = Y_Ca44_Sc44;
plot(x,Y_C11b,'k'); hold on
plot(x,Y_N13b,'c')
plot(x,Y_O15b,'m')
plot(x,Y_C10b,'y');
plot(x,Y_Sc44b,'g');
legend('C11','N13','O15','C10','Sc44','Location', 'northwest');
[f,g]=min(Ddepb);
axis([0 (ceil(g*dx)) 0 max(max(Y_C11b,Y_O15b)+0.2*max(Y_C11b,Y_O15b))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

% Figure in PMMA
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('PMMA');
hold on
plot(x,100*Ddepp)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepp)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11p = Y_O16_C11p + Y_C12_C11p + Y_N14_C11p;
Y_N13p = Y_O16_N13p + Y_N14_N13p;
Y_O15p = Y_O16_O15p;
Y_C10p = Y_C12_C10p;
hold on
plot(x,Y_C11p,'k');
plot(x,Y_O15p,'m'); 
plot(x,Y_N13p,'c');
plot(x,Y_C10p,'y');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('C11','O15','N13','C10','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 (max(Y_C11p+Y_O15p)+0.2*max(Y_C11p+Y_O15p))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

% Figure in PMMA
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('PMMA');
hold on
plot(x,100*Ddepp)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepp)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11p = Y_O16_C11p + Y_C12_C11p + Y_N14_C11p;
Y_N13p = Y_O16_N13p + Y_N14_N13p;
Y_O15p = Y_O16_O15p;
Y_C10p = Y_C12_C10p;
plot(x,(Y_C11p+Y_O15p),'b'); hold on
plot(x,Y_C11p,'k');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('Total','C11','N13','O15','C10','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 (max(Y_C11p+Y_O15p)+0.2*max(Y_C11p+Y_O15p))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

%% PG

%PG PMMA
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('PMMA');
hold on
plot(x,100*Ddepp)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepp)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4p = Y_PG_C12_C12_4p+Y_PG_O16_C12_4p;
Y_1p = Y_PG_N14_N14_1p;
Y_6p = Y_PG_O16_O16_6p;
plot(x,Y_4p,'b'); hold on
plot(x,Y_6p,'m');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('4.44 MeV', '6.13 MeV','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 (max(Y_4p)+0.2*max(Y_4p))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

%PG Piel
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Tissue');
hold on
plot(x,100*Ddept)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddept)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4t = Y_PG_C12_C12_4t+Y_PG_O16_C12_4t;
Y_1t = Y_PG_N14_N14_1t;
Y_2t = Y_PG_N14_N14_2t;
Y_6t = Y_PG_O16_O16_6t;
plot(x,Y_6t,'b'); hold on
plot(x,Y_4t,'m');
plot(x,Y_2t,'c');
plot(x,Y_1t,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddept);
axis([0 (ceil(g*dx)) 0 (max(Y_4t)+0.2*max(Y_4t))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

%PG Adipose
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Adipose');
hold on
plot(x,100*Ddepa)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepa)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4a = Y_PG_C12_C12_4a+Y_PG_O16_C12_4a;
Y_1a = Y_PG_N14_N14_1a;
Y_2a = Y_PG_N14_N14_2a;
Y_6a = Y_PG_O16_O16_6a;
plot(x,Y_6a,'b'); hold on
plot(x,Y_4a,'m');
plot(x,Y_2a,'c');
plot(x,Y_1a,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepa);
axis([0 (ceil(g*dx)) 0 (max(Y_4a)+0.2*max(Y_4a))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

%PG Bone
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Bone');
hold on
plot(x,100*Ddepb)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepb)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4b = Y_PG_C12_C12_4b+Y_PG_O16_C12_4b;
Y_1b = Y_PG_N14_N14_1b;
Y_2b = Y_PG_N14_N14_2b;
Y_6b = Y_PG_O16_O16_6b;
plot(x,Y_6b,'b'); hold on
plot(x,Y_4b,'m');
plot(x,Y_2b,'c');
plot(x,Y_1b,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepb);
axis([0 (ceil(g*dx)) 0 (max(Y_4b)+0.2*max(Y_4b))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)


%% PET activity (t)
calcTimes = [0 60 1000 3600]; % s
act_C11 = zeros(numel(calcTimes), numel(x));
act_N13 = zeros(numel(calcTimes), numel(x));
act_O15 = zeros(numel(calcTimes), numel(x));

act_C11t = zeros(numel(calcTimes), numel(x));
act_C10t = zeros(numel(calcTimes), numel(x));
act_N13t = zeros(numel(calcTimes), numel(x));
act_O15t= zeros(numel(calcTimes), numel(x));

act_C11a = zeros(numel(calcTimes), numel(x));
act_C10a = zeros(numel(calcTimes), numel(x));
act_N13a = zeros(numel(calcTimes), numel(x));
act_O15a= zeros(numel(calcTimes), numel(x));

act_C11b = zeros(numel(calcTimes), numel(x));
act_C10b = zeros(numel(calcTimes), numel(x));
act_N13b = zeros(numel(calcTimes), numel(x));
act_O15b = zeros(numel(calcTimes), numel(x));
act_Sc44b = zeros(numel(calcTimes), numel(x));

act_C11p = zeros(numel(calcTimes), numel(x));
act_C10p = zeros(numel(calcTimes), numel(x));
act_N13p = zeros(numel(calcTimes), numel(x));
act_O15p = zeros(numel(calcTimes), numel(x));

for i=1:numel(calcTimes)
    % Water
    act_C11(i,:) = landa_C11 .* Y_O16_C11s .* exp(- landa_C11 * calcTimes(i));
    act_N13(i,:) = landa_N13 .* Y_O16_N13s .* exp(- landa_N13 * calcTimes(i));
    act_O15(i,:) = landa_O15 .* Y_O16_O15s .* exp(- landa_O15 * calcTimes(i));
    
    % Tissue
    act_C11t(i,:) = landa_C11 .* Y_C11t .* exp(- landa_C11 * calcTimes(i));
    act_C10t(i,:) = landa_C10 .* Y_C10t .* exp(- landa_C10 * calcTimes(i));    
    act_N13t(i,:) = landa_N13 .* Y_N13t .* exp(- landa_N13 * calcTimes(i));
    act_O15t(i,:) = landa_O15 .* Y_O15t .* exp(- landa_O15 * calcTimes(i));  
    
     % Adipose
    act_C11a(i,:) = landa_C11 .* Y_C11a .* exp(- landa_C11 * calcTimes(i));
    act_C10a(i,:) = landa_C10 .* Y_C10a .* exp(- landa_C10 * calcTimes(i));    
    act_N13a(i,:) = landa_N13 .* Y_N13a .* exp(- landa_N13 * calcTimes(i));
    act_O15a(i,:) = landa_O15 .* Y_O15a .* exp(- landa_O15 * calcTimes(i)); 
    
     % Bone
    act_C11b(i,:) = landa_C11 .* Y_C11b .* exp(- landa_C11 * calcTimes(i));
    act_C10b(i,:) = landa_C10 .* Y_C10b .* exp(- landa_C10 * calcTimes(i));    
    act_N13b(i,:) = landa_N13 .* Y_N13b .* exp(- landa_N13 * calcTimes(i));
    act_O15b(i,:) = landa_O15 .* Y_O15b .* exp(- landa_O15 * calcTimes(i)); 
    act_Sc44b(i,:) = landa_Sc44 .* Y_Sc44b .* exp(- landa_Sc44 * calcTimes(i));
    
     % Bone
    act_C11p(i,:) = landa_C11 .* Y_C11p .* exp(- landa_C11 * calcTimes(i));
    act_C10p(i,:) = landa_C10 .* Y_C10p .* exp(- landa_C10 * calcTimes(i));    
    act_N13p(i,:) = landa_N13 .* Y_N13p .* exp(- landa_N13 * calcTimes(i));
    act_O15p(i,:) = landa_O15 .* Y_O15p .* exp(- landa_O15 * calcTimes(i)); 
end
act_total = (act_C11 + act_N13 + act_O15);
act_totalt = (act_C11t + act_C10t + act_N13t + act_O15t);
act_totala = (act_C11a + act_C10a + act_N13a + act_O15a);
act_totalb = (act_C11b + act_C10b + act_N13b + act_O15b+act_Sc44b);
act_totalp = (act_C11p + act_C10p + act_N13p + act_O15p);

F=figure;
G=figure;
H=figure;
I=figure;
AAA=length(calcTimes);

for i=1:numel(calcTimes)
    figure(F)
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totalt(i,:),'b')
    hold on
    plot(x, act_C11t(i,:),'r')
    plot(x, act_N13t(i,:),'y')
    plot(x, act_O15t(i,:),'c')
    plot(x, act_C10t(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,100*Ddept);
    ylabel('-dE/dx')
    [f,g]=min(Ddept);
    axis([0 (ceil(g*dx)) 0 max(100*Ddept)+0.2*max(100*Ddept)]);
    
    figure(G)
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totala(i,:),'b')
    hold on
    plot(x, act_C11a(i,:),'r')
    plot(x, act_N13a(i,:),'y')
    plot(x, act_O15a(i,:),'c')
    plot(x, act_C10a(i,:),'g')
    title(sprintf('Activity at t=%i s',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
    axis([0 5 0 max(act_totala(i,:))+0.2*max(act_totala(i,:))]);   
    yyaxis right
    grid on
    plot(x,100*Ddepa);
    ylabel('-dE/dx')
    [f,g]=min(Ddepa);
    axis([0 (ceil(g*dx)) 0 max(100*Ddepa)+0.2*max(100*Ddepa)]);
    
    figure(H)
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totalb(i,:),'b')
    hold on
    plot(x, act_C11b(i,:),'r')
    plot(x, act_N13b(i,:),'y')
    plot(x, act_O15b(i,:),'c')
    plot(x, act_C10b(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16) 
    axis([0 11 0 max(act_totalb(i,:))+0.2*max(act_totalb(i,:))]);     
    yyaxis right
    grid on
    plot(x,100*Ddepb);
    ylabel('-dE/dx')
    [f,g]=min(Ddepb);
    axis([0 (ceil(g*dx)) 0 max(100*Ddepb)+0.2*max(100*Ddepb)]);
    
    figure(I);
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totalp(i,:),'b')
    hold on
    plot(x, act_C11p(i,:),'r')
    plot(x, act_N13p(i,:),'y')
    plot(x, act_O15p(i,:),'c')
    plot(x, act_C10p(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,100*Ddepp);
    ylabel('-dE/dx')
    [f,g]=min(Ddepp);
    axis([0 (ceil(g*dx)) 0 max(100*Ddepp)+0.2*max(100*Ddepp)]);
end



%% Calculo a Dosis

%Ddept=Ddept/rho_tissue;
%Ddepa=Ddepa/rho_adipose;
%Ddepb=Ddepb/rho_bone;
%Ddepp=Ddepp/rho_PMMA;
%Dosist=0;
%Dosisa=0;
%Dosisb=0;

%for i=1:(numel(x)-1)
    
    %Dosist=Dosist+dx*(Ddept(i)+Ddept(i+1))/2; %MeV cm2/g
    %Dosisa=Dosisa+dx*(Ddepa(i)+Ddepa(i+1))/2; %MeV cm2/g
    %Dosisb=Dosisb+dx*(Ddepb(i)+Ddepb(i+1))/2; %MeV cm2/g
    



%end


%Dosist=Dosist*1.6e-7; %Gy cm2
%Dosisa=Dosisa*1.6e-7; %Gy cm2
%Dosisb=Dosisb*1.6e-7; %Gy cm2

%Dosist=max(Ddept)*1.6e-7; %Gy cm2
%Dosisa=max(Ddepa)*1.6e-7; %Gy cm2
%Dosisb=max(Ddepb)*1.6e-7; %Gy cm2
%Numt=1/Dosist;
%Numa=1/Dosisa;
%Numb=1/Dosisb;


%% Actividad con el tiempo
%Por ahora suponemos un protón por segundo.
deltat=1;
a=120;
c=1;
t=900;
temp_C10t=zeros(t+1,numel(x));temp_C10a=zeros(t+1,numel(x));
temp_C11t=zeros(t+1,numel(x));temp_C11a=zeros(t+1,numel(x));
temp_N13t=zeros(t+1,numel(x));temp_N13a=zeros(t+1,numel(x));
temp_O15t=zeros(t+1,numel(x));temp_O15a=zeros(t+1,numel(x));
temp_C10p=zeros(t+1,numel(x));temp_C10b=zeros(t+1,numel(x));
temp_C11p=zeros(t+1,numel(x));temp_C11b=zeros(t+1,numel(x));
temp_N13p=zeros(t+1,numel(x));temp_N13b=zeros(t+1,numel(x));
temp_O15p=zeros(t+1,numel(x));temp_O15b=zeros(t+1,numel(x));
temp_totalp2=zeros(t+1);
temp_totalt2=zeros(t+1);
temp_totalb2=zeros(t+1);
temp_totala2=zeros(t+1);
temp_parcC11t=zeros(t+1);temp_parcC11p=zeros(t+1);temp_parcC11a=zeros(t+1);temp_parcC11b=zeros(t+1);
temp_parcO15t=zeros(t+1);temp_parcO15p=zeros(t+1);temp_parcO15a=zeros(t+1);temp_parcO15b=zeros(t+1);
temp_parcN13t=zeros(t+1);temp_parcN13p=zeros(t+1);temp_parcN13a=zeros(t+1);temp_parcN13b=zeros(t+1);
temp_parcC10t=zeros(t+1);temp_parcC10p=zeros(t+1);temp_parcC10a=zeros(t+1);temp_parcC10b=zeros(t+1);
for i=1:t+1;
    b=i;
    if i<a
        for j=1:b
            
            temp_C10t(i,:)=temp_C10t(i,:)+deltat * Y_C10t.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11t(i,:)=temp_C11t(i,:)+deltat * Y_C11t.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13t(i,:)=temp_N13t(i,:)+deltat * Y_N13t.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15t(i,:)=temp_O15t(i,:)+deltat * Y_O15t.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
            
            temp_C10p(i,:)=temp_C10p(i,:)+deltat * Y_C10p.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11p(i,:)=temp_C11p(i,:)+deltat * Y_C11p.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13p(i,:)=temp_N13p(i,:)+deltat * Y_N13p.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15p(i,:)=temp_O15p(i,:)+deltat * Y_O15p.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
         
            temp_C10a(i,:)=temp_C10a(i,:)+deltat * Y_C10a.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11a(i,:)=temp_C11a(i,:)+deltat * Y_C11a.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13a(i,:)=temp_N13a(i,:)+deltat * Y_N13a.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15a(i,:)=temp_O15a(i,:)+deltat * Y_O15a.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
            
            temp_C10b(i,:)=temp_C10b(i,:)+deltat * Y_C10b.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11b(i,:)=temp_C11b(i,:)+deltat * Y_C11b.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13b(i,:)=temp_N13b(i,:)+deltat * Y_N13b.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15b(i,:)=temp_O15b(i,:)+deltat * Y_O15b.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
            
        end
    else
            for j=1:a;
            temp_C10t(i,:)=temp_C10t(i,:)+deltat *Y_C10t.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11t(i,:)=temp_C11t(i,:)+deltat *Y_C11t.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13t(i,:)=temp_N13t(i,:)+deltat *Y_N13t.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15t(i,:)=temp_O15t(i,:)+deltat *Y_O15t.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c));
            
            temp_C10p(i,:)=temp_C10p(i,:)+deltat *Y_C10p.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11p(i,:)=temp_C11p(i,:)+deltat *Y_C11p.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13p(i,:)=temp_N13p(i,:)+deltat *Y_N13p.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15p(i,:)=temp_O15p(i,:)+deltat *Y_O15p.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c));
            
            temp_C10a(i,:)=temp_C10a(i,:)+deltat *Y_C10a.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11a(i,:)=temp_C11a(i,:)+deltat *Y_C11a.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13a(i,:)=temp_N13a(i,:)+deltat *Y_N13a.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15a(i,:)=temp_O15a(i,:)+deltat *Y_O15a.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c));
            
            temp_C10b(i,:)=temp_C10b(i,:)+deltat *Y_C10b.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11b(i,:)=temp_C11b(i,:)+deltat *Y_C11b.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13b(i,:)=temp_N13b(i,:)+deltat *Y_N13b.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15b(i,:)=temp_O15b(i,:)+deltat *Y_O15b.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c)); 
            
            end
            c=c+1;
    end
           
            
   % else;
           % temp_C10t(i,:)=temp_C10t(i,:)+Y_C12_C10t.*landa_C10.*exp(-landa_C10*(t-1));
            
            
            
     
    temp_totalt2(i)=sum(temp_C10t(i,:))+sum(temp_C11t(i,:))+sum(temp_N13t(i,:))+sum(temp_O15t(i,:));
    temp_parcC11t(i)=sum(temp_C11t(i,:));
    temp_parcN13t(i)=sum(temp_N13t(i,:));
    temp_parcC10t(i)=sum(temp_C10t(i,:));
    temp_parcO15t(i)=sum(temp_O15t(i,:));
    
    temp_totalp2(i)=sum(temp_C10p(i,:))+sum(temp_C11p(i,:))+sum(temp_N13p(i,:))+sum(temp_O15p(i,:));
    temp_parcC11p(i)=sum(temp_C11p(i,:));
    temp_parcN13p(i)=sum(temp_N13p(i,:));
    temp_parcC10p(i)=sum(temp_C10p(i,:));
    temp_parcO15p(i)=sum(temp_O15p(i,:));
    
    temp_totala2(i)=sum(temp_C10a(i,:))+sum(temp_C11a(i,:))+sum(temp_N13a(i,:))+sum(temp_O15a(i,:));
    temp_parcC11a(i)=sum(temp_C11a(i,:));
    temp_parcN13a(i)=sum(temp_N13a(i,:));
    temp_parcC10a(i)=sum(temp_C10a(i,:));
    temp_parcO15a(i)=sum(temp_O15a(i,:));
    
    temp_totalb2(i)=sum(temp_C10b(i,:))+sum(temp_C11b(i,:))+sum(temp_N13b(i,:))+sum(temp_O15b(i,:));
    temp_parcC11b(i)=sum(temp_C11b(i,:));
    temp_parcN13b(i)=sum(temp_N13b(i,:));
    temp_parcC10b(i)=sum(temp_C10b(i,:));
    temp_parcO15b(i)=sum(temp_O15b(i,:));

    
    
end
    temp_totalt=temp_C10t+temp_C11t+temp_N13t+temp_O15t;
    temp_totalp=temp_C10p+temp_C11t+temp_N13p+temp_O15p;
    temp_totala=temp_C10a+temp_C11a+temp_N13a+temp_O15a;
    temp_totalb=temp_C10b+temp_C11b+temp_N13b+temp_O15b;
    




    %% Dibuja la Actividad total en función del tiempo.
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totalt2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15t(:,1));
    plot(T*deltat,temp_parcC11t(:,1));
    plot(T*deltat,temp_parcN13t(:,1));
    plot(T*deltat,temp_parcC10t(:,1));
    title('Actividad total en a lo largo del tiempo HUESO');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalt2(:,1))+0.2*max(temp_totalt2(:,1)))]); 
    
    figure;
    plot(T*deltat,(24*temp_totalp2(:,1)));
    hold on;
    plot(T*deltat,24*temp_parcO15p(:,1));
    plot(T*deltat,24*temp_parcC11p(:,1));
    plot(T*deltat,24*temp_parcN13p(:,1));
    plot(T*deltat,24*temp_parcC10p(:,1));
%    title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(24*temp_totalp2(:,1))+0.2*max(24*temp_totalp2(:,1)))]);
    
    figure;
    plot(T*deltat,(temp_totala2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15a(:,1));
    plot(T*deltat,temp_parcC11a(:,1));
    plot(T*deltat,temp_parcN13a(:,1));
    plot(T*deltat,temp_parcC10a(:,1));
    title('Actividad total en a lo largo del tiempo ADIPOSE');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(temp_totala2(:,1))+0.2*max(temp_totala2(:,1)))]);
    
    figure;
    plot(T*deltat,(temp_totalb2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15b(:,1));
    plot(T*deltat,temp_parcC11b(:,1));
    plot(T*deltat,temp_parcN13b(:,1));
    plot(T*deltat,temp_parcC10b(:,1));
    title('Actividad total en a lo largo del tiempo HUESO');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(temp_totalb2(:,1))+0.2*max(temp_totalb2(:,1)))]);

%% Actividad total en función de tiempo
%Por ahora suponemos un protón por segundo.
cont1=0;
cont2=0;
c=1;
int_totalp2=zeros(t+1);
int_totalt2=zeros(t+1);
int_parcC11t=zeros(t+1);
int_parcO15t=zeros(t+1);
int_parcN13t=zeros(t+1);
int_parcC10t=zeros(t+1);
int_parcC11p=zeros(t+1);
int_parcO15p=zeros(t+1);
int_parcN13p=zeros(t+1);
int_parcC10p=zeros(t+1);
int_C11t=zeros(t+1,numel(x));
int_O15t=zeros(t+1,numel(x));
int_N13t=zeros(t+1,numel(x));
int_C10t=zeros(t+1,numel(x));
int_C11p=zeros(t+1,numel(x));
int_O15p=zeros(t+1,numel(x));
int_N13p=zeros(t+1,numel(x));
int_C10p=zeros(t+1,numel(x));
for i=2:t+1;
    b=i;
    if i<a
        cont1=cont1+1;
        for j=1:b;
            
            %Tissue
            int_C10t(i,:)=int_C10t(i,:)+deltat*Y_C10t.*(-exp(-landa_C10*(deltat*j))-exp(-landa_C10*(deltat*j-deltat*1)));
            int_C11t(i,:)=int_C11t(i,:)+deltat*Y_C11t.*(-exp(-landa_C11*(deltat*j))+exp(-landa_C11*(deltat*j-deltat*1)));
            int_N13t(i,:)=int_N13t(i,:)+deltat*Y_N13t.*(-exp(-landa_N13*(deltat*j))+exp(-landa_N13*(deltat*j-deltat*1)));;
            int_O15t(i,:)=int_O15t(i,:)+deltat*Y_O15t.*(-exp(-landa_O15*(deltat*j))+exp(-landa_O15*(deltat*j-deltat*1)));;
            
            %PMMA
            int_C10p(i,:)=int_C10p(i,:)+deltat*Y_C10p.*(-exp(-landa_C10*(deltat*j))+exp(-landa_C10*(deltat*j-deltat*1)));
            int_C11p(i,:)=int_C11p(i,:)+deltat*Y_C11p.*(-exp(-landa_C11*(deltat*j))+exp(-landa_C11*(deltat*j-deltat*1)));
            int_N13p(i,:)=int_N13p(i,:)+deltat*Y_N13p.*(-exp(-landa_N13*(deltat*j))+exp(-landa_N13*(deltat*j-deltat*1)));
            int_O15p(i,:)=int_O15p(i,:)+deltat*Y_O15p.*(-exp(-landa_O15*(deltat*j))+exp(-landa_O15*(deltat*j-deltat*1)));
           
        end

    else
       
           cont2=cont2+1;
        for j=1:a;
            
            %Tissue
            int_C10t(i,:)=int_C10t(i,:)+deltat*Y_C10t.*(-exp(-landa_C10*(deltat*j+deltat*c))+exp(-landa_C10*(deltat*j-deltat*1+deltat*c)));
            int_C11t(i,:)=int_C11t(i,:)+deltat*Y_C11t.*(-exp(-landa_C11*(deltat*j+deltat*c))+exp(-landa_C11*(deltat*j-deltat*1+deltat*c)));
            int_N13t(i,:)=int_N13t(i,:)+deltat*Y_N13t.*(-exp(-landa_N13*(deltat*j+deltat*c))+exp(-landa_N13*(deltat*j-deltat*1+deltat*c)));
            int_O15t(i,:)=int_O15t(i,:)+deltat*Y_O15t.*(-exp(-landa_O15*(deltat*j+deltat*c))+exp(-landa_O15*(deltat*j-deltat*1+deltat*c)));
            
            %PMMA
            int_C10p(i,:)=int_C10p(i,:)+deltat*Y_C10p.*(-exp(-landa_C10*(deltat*j+deltat*c))+exp(-landa_C10*(deltat*j-deltat*1+deltat*c)));
            int_C11p(i,:)=int_C11p(i,:)+deltat*Y_C11p.*(-exp(-landa_C11*(deltat*j+deltat*c))+exp(-landa_C11*(deltat*j-deltat*1+deltat*c)));
            int_N13p(i,:)=int_N13p(i,:)+deltat*Y_N13p.*(-exp(-landa_N13*(deltat*j+deltat*c))+exp(-landa_N13*(deltat*j-deltat*1+deltat*c)));
            int_O15p(i,:)=int_O15p(i,:)+deltat*Y_O15p.*(-exp(-landa_O15*(deltat*j+deltat*c))+exp(-landa_O15*(deltat*j-deltat*1+deltat*c)));
           
        end
            c=c+1;
            

      
    end
           
    int_totalt2(i)=sum(int_C10t(i,:))+sum(int_C11t(i,:))+sum(int_N13t(i,:))+sum(int_O15t(i,:));
    int_parcC11t(i)=sum(int_C11t(i,:));
    int_parcO15t(i)=sum(int_O15t(i,:));
    int_parcN13t(i)=sum(int_N13t(i,:));
    int_parcC10t(i)=sum(int_C10t(i,:));
    
    int_totalp2(i)=sum(int_C10p(i,:))+sum(int_C11p(i,:))+sum(int_N13p(i,:))+sum(int_O15p(i,:));
    int_parcC11p(i)=sum(int_C11p(i,:));
    int_parcO15p(i)=sum(int_O15p(i,:));
    int_parcN13p(i)=sum(int_N13p(i,:));
    int_parcC10p(i)=sum(int_C10p(i,:));

    
end

%% Actividadd total 2

    figure;
    plot(T*deltat,(24.*int_totalp2(:,1)));
    hold on
    plot(T*deltat,24.*int_parcO15p(:,1))
    plot(T*deltat,24.*int_parcC11p(:,1))
    plot(T*deltat,24.*int_parcC10p(:,1))
    plot(T*deltat,24.*int_parcN13p(:,1))
    title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)')
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','C10','N13','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    axis([0 t*deltat 0 (max(24.*int_totalp2(:,1))+0.2*max(24.*int_totalp2(:,1)))]);

    
    
%% Actividad durante Bean-On

beam=zeros(1,length(x));
beam_onp=zeros(1,length(x));
beam_offp=zeros(1,length(x));
C10_onp=zeros(1,length(x));
C10_offp=zeros(1,length(x));
C11_onp=zeros(1,length(x));
C11_offp=zeros(1,length(x));
N13_onp=zeros(1,length(x));
N13_offp=zeros(1,length(x));
O15_onp=zeros(1,length(x));
O15_offp=zeros(1,length(x));

for i=1:240
    beam=beam+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
    if i<a
        beam_onp=beam_onp+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
        C10_onp=C10_onp+int_C10p(i,:);
        C11_onp=C11_onp+int_C11p(i,:);
        N13_onp=N13_onp+int_N13p(i,:);
        O15_onp=O15_onp+int_O15p(i,:);
    else
        beam_offp=beam_offp+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
        C10_offp=C10_offp+int_C10p(i,:);
        C11_offp=C11_offp+int_C11p(i,:);
        N13_offp=N13_offp+int_N13p(i,:);
        O15_offp=O15_offp+int_O15p(i,:);
    end
end

figure;

%subplot(2,2,1);
    yyaxis right
    grid on
    plot(x,100*Ddepp,'k--');
    ylabel('-dE/dx')
    axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
    yyaxis left
hold on;
plot(x,24*beam_onp(1,:),'b');
plot(x,24*beam_offp(1,:),'r');
    title('Actividad total ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Beam On','Beam Off','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (24*max(beam_offp)+0.2*max(24*beam_offp))]);

figure;
%subplot(2,2,2);
    yyaxis right
    grid on
    plot(x,100*Ddepp,'k--');
    ylabel('-dE/dx')
    axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
    yyaxis left
hold on;
plot(x,C11_onp(1,:),'b');
plot(x,C11_offp(1,:),'r');
    title('C11 ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Beam On','Beam Off','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(C11_offp)+0.2*max(C11_offp))]);

figure;
%subplot(2,2,3);
    yyaxis right
    grid on
    plot(x,100*Ddepp,'k--');
    ylabel('-dE/dx')
    axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
    yyaxis left
hold on;
plot(x,O15_onp(1,:),'b');
plot(x,O15_offp(1,:),'r');
    title('O15 ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Beam On','Beam Off','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(O15_offp)+0.2*max(O15_offp))]);

figure;
%subplot(2,2,4);
    yyaxis right
    grid on
    plot(x,100*Ddepp,'k--');
    ylabel('-dE/dx')
    axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
    yyaxis left
hold on;
plot(x,N13_onp(1,:),'b');
plot(x,N13_offp(1,:),'r');
    title('N13 ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Beam On','Beam Off','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(N13_offp)+0.2*max(N13_offp))]);

figure
    yyaxis right
    grid on
    plot(x,100*Ddepp,'k--');
    ylabel('-dE/dx')
    axis([0 5 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
    yyaxis left
hold on;
plot(x,beam(1,:),'b');
    title('Total ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Beam On','Dose','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(beam(1,:))+0.2*max(beam(1,:)))]);



