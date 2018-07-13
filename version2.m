% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es

%% Par?metros a modificar:
clear all;close all;

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control1.mat');
load('control2.mat');

%PARAMETROS
dx=0.1        %Paso del intervalo 
xref=20;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=60;        %Energía inicial del haz
deltat=1;
a=60/deltat;
t=900/deltat;

pps=1;         %protones/segundo
%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
PMMA_Molar=100.12; %g/mol
rho_w = 1; % g/cm3
rho_Zn = 7.14; %g/cm3
rho_tissue = 1.1; %g/cm3
rho_bone = 1.85; %g/cm3
rho_adipose = 0.92; %g/cm3
rho_PMMA= 1.18; %g/cm3
rho_cartilage = 1.1; %g/cm3
rho_muscle = 1.04; %g/cm3
%% 

%Densidades Atomicas (A falta de las masas molares usamos
%sum(Comp_bone.*W_ele)) CORREGIR
%rho_tissue_A = 4.6243E+22; % atoms/cm3
rho_tissue_A = AvNmbr*rho_tissue/sum(Comp_tissue.*W_ele);  % atoms/cm3
rho_bone_A = AvNmbr*rho_bone/sum(Comp_bone.*W_ele);  % atoms/cm3
rho_adipose_A = AvNmbr*rho_adipose/sum(Comp_adipose.*W_ele);  % atoms/cm3
rho_cartilage_A = AvNmbr*rho_cartilage/sum(Comp_cartilage.*W_ele); %atoms/cm3
rho_muscle_A = AvNmbr*rho_muscle/sum(Comp_muscle.*W_ele); % atoms/cm3
%rho_PMMA_A = AvNmbr*rho_PMMA/sum(Comp_PMMA.*W_ele);  % atoms/cm3
rho_PMMA_A = AvNmbr*rho_PMMA/PMMA_Molar;  % atoms/cm3
rho_w_A =  rho_w * AvNmbr / waterMolecularWeight; % molecules / cm3
ZnAtomicWeight = 65.38; % g/mol
rho_Zn_A = Zn_fraction * rho_Zn * AvNmbr / ZnAtomicWeight; % molecules / cm3

%% 

%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water () muscle(m) cartilage(c)
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
rho_O16_Am = rho_muscle_A * Comp_muscle(4) * O16_ab;
rho_N14_Am = rho_muscle_A * Comp_muscle(3) * N14_ab;
rho_C12_Am = rho_muscle_A * Comp_muscle(2) * C12_ab;
rho_O16_Ac = rho_cartilage_A * Comp_cartilage(4) * O16_ab;
rho_N14_Ac = rho_cartilage_A * Comp_cartilage(3) * N14_ab;
rho_C12_Ac = rho_cartilage_A * Comp_cartilage(2) * C12_ab;

%% 

%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
E = nan(size(x));
Et = nan(size(x));
Eb = nan(size(x));
Ea = nan(size(x));
Ep = nan(size(x));
Ec = nan(size(x));
Em = nan(size(x));

%Energia depositada por material
Ddep = zeros(size(E));
Ddept = zeros(size(E));
Ddepb = zeros(size(E));
Ddepa = zeros(size(E));
Ddepp = zeros(size(E));
Ddepc = zeros(size(E));
Ddepm = zeros(size(E));

%Energia actual de cada material
currentE = E0;
currentEt = E0;
currentEb = E0;
currentEa = E0;
currentEp = E0;
currentEc = E0;
currentEm = E0;

%CREACION DE VECTORES DE YIELD
% In water (full + simplified versions)
Y_O16_C11 = nan(size(x));
Y_O16_N13 = nan(size(x));
Y_O16_O15 = nan(size(x));
Y_O16_C11s = nan(size(x));
Y_O16_N13s = nan(size(x));
Y_O16_O15s = nan(size(x));
Y_PG_C12_C12_4w = zeros(size(x));
Y_PG_O16_C12_4w = zeros(size(x));
Y_PG_N14_N14_1w = zeros(size(x));
Y_PG_N14_N14_2w = zeros(size(x));
Y_PG_O16_O16_6w = zeros(size(x));

% In tissue (simplified form only)
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
Y_PG_trap = zeros(size(x));

% In cartilage (simplified form only)
Y_O16_C11c = zeros(size(x));
Y_O16_N13c = zeros(size(x));
Y_O16_O15c = zeros(size(x));
Y_N14_N13c = zeros(size(x));
Y_N14_O14c = zeros(size(x));
Y_N14_C11c = zeros(size(x));   
Y_C12_C11c = zeros(size(x));     
Y_C12_C10c = zeros(size(x)); 
Y_PG_C12_C12_4c = zeros(size(x));
Y_PG_O16_C12_4c = zeros(size(x));
Y_PG_N14_N14_1c = zeros(size(x));
Y_PG_N14_N14_2c = zeros(size(x));
Y_PG_O16_O16_6c = zeros(size(x));

% In muscle (simplified form only)
Y_O16_C11m = zeros(size(x));
Y_O16_N13m = zeros(size(x));
Y_O16_O15m = zeros(size(x));
Y_N14_N13m = zeros(size(x));
Y_N14_O14m = zeros(size(x));
Y_N14_C11m = zeros(size(x));   
Y_C12_C11m = zeros(size(x));     
Y_C12_C10m = zeros(size(x)); 
Y_PG_C12_C12_4m = zeros(size(x));
Y_PG_O16_C12_4m = zeros(size(x));
Y_PG_N14_N14_1m = zeros(size(x));
Y_PG_N14_N14_2m = zeros(size(x));
Y_PG_O16_O16_6m = zeros(size(x));

%% 

%CALCULO YIELD
%Recorremos cada una de las regiones del espacio caculando analíticamente
%el número de protones generado en cada una de ellas.


for i=1:(numel(x)-1)
    
    %Calculamos el poder de frenado para la energía con la que llega el
    %protón a la lámina.
    S_w = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Zn = max(0,1000*S_Zn_F(currentE*1000)); % MeV/(g/cm2) 
    S_Znt = max(0,1000*S_Zn_F(currentEt*1000)); % MeV/(g/cm2) 
    S_t = max(0,1000*S_t_F(currentEt*1000)); % MeV/(g/cm2)
    S_b = max(0,1000*S_b_F(currentEb*1000)); 
    S_a = max(0,1000*S_a_F(currentEa*1000));
    S_p = max(0,1000*S_p_F(currentEp*1000));
    S_c = max(0,1000*S_c_F(currentEc*1000));
    S_m = max(0,1000*S_m_F(currentEm*1000));
    
    %Multiplicamospor la densidad al poder de frenado
    S1 = (S_w*rho_w); % MeV/cm
    S1t = (S_t*rho_tissue); %*(1-Zn_fraction) + S_Znt*rho_Zn*Zn_fraction); % MeV/cm
    S1b = (S_b*rho_bone); % MeV/cm
    S1a = (S_a*rho_adipose); % MeV/cm
    S1p = (S_p*rho_PMMA); % MeV/cm
    S1c = (S_c*rho_cartilage); % MeV/cm
    S1m = (S_m*rho_muscle); %MeV/cm   
    %Calculamos la energía que pierde el protón al atravesar la lámina (se
    %supone que toda la pérdida se hace al final).
    deltaE = dx*S1; % MeV
    deltaEt = dx*S1t; % MeV
    deltaEb = dx*S1b; % MeV
    deltaEa = dx*S1a; %MeV
    deltaEp = dx*S1p; %MeV
    deltaEc = dx*S1c; %MeV
    deltaEm = dx*S1m; %MeV
    
    %Guardamos las energías para poder representar todo luego en función de
    %ella.
    E(i) = currentE; % MeV
    Et(i) = currentEt; % MeV
    Eb(i) = currentEb;  % MeV
    Ea(i) = currentEa; % MeV
    Ep(i) = currentEp; % MeV
    Ec(i) = currentEc; %MeV
    Em(i) = currentEm; %MeV
    
    %Calculamos la energía con la que sale de la lámina que se usará en la
    %siguiente iteración.
    currentE = currentE - deltaE; % MeV
    currentEt = currentEt - deltaEt; % MeV
    currentEb = currentEb - deltaEb; % MeV
    currentEa = currentEa - deltaEa; %MeV
    currentEp = currentEp - deltaEp; %MeV
    currentEc = currentEc - deltaEc; %MeV
    currentEm = currentEm - deltaEm; %MeV
    
    %NO SE USA
    S_w2 = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Zn2 = max(0,1000*S_Zn_F(currentE*1000)); % MeV/(g/cm2)
    S_t2 = max(0,1000*S_t_F(currentEt*1000)); % MeV/(g/cm2)
    S_Zn2t = max(0,1000*S_Zn_F(currentEt*1000)); % MeV/(g/cm2)
    S_b2 = max(0,1000*S_b_F(currentEb*1000));% MeV/(g/cm2)
    S_a2 = max(0,1000*S_a_F(currentEa*1000));% MeV/(g/cm2)
    S_p2 = max(0,1000*S_p_F(currentEp*1000));% MeV/(g/cm2)
    S_c2 = max(0,1000*S_c_F(currentEc*1000)); %MeV/(g/cm2)
    S_m2 = max(0,1000*S_c_F(currentEm*1000));
    S2 = (S_w2*rho_w); % MeV/cm
    S2t = (S_t2*rho_tissue); %*(1-Zn_fraction) + S_Zn2t*rho_Zn*Zn_fraction); % MeV/cm
    S2b = (S_b2*rho_bone); % MeV/cm3
    S2a = (S_a2*rho_adipose); % MeV/cm3
    S2p = (S_p2*rho_PMMA); % MeV/cm3
    S2c = (S_c2*rho_cartilage); %MeV/cm3 
    S2m = (S_m2*rho_muscle);
    
    %Guardamos la energía depositada en cada uno de los materiales.
    Ddep(i) = deltaE; % MeV
    Ddept(i) = deltaEt; % MeV
    Ddepb(i) = deltaEb;  % MeV
    Ddepa(i) = deltaEa; %MeV
    Ddepp(i) = deltaEp; %MeV
    Ddepc(i) = deltaEc; %MeV
    Ddepm(i) = deltaEm; %MeV
    
    %Creamos unas variables con los datos anteriores para introducir a
    %continuación en los Yield.
    E1 = E(i);
    E1t = Et(i);
    E1b = Eb(i);
    E1a = Ea(i);
    E1p = Ep(i);
    E1c = Ec(i);
    E1m = Em(i);
    E2 = currentE;
    E2t = currentEt;
    E2a = currentEa;
    E2b = currentEb;
    E2p = currentEp;
    E2c = currentEc;
    E2m = currentEm;
    
    %Estas de abajo solo se usan por el método trapezoidal que no es el que
    %estamos teniendo en cuenta.
    E2E1 = [currentE E(i)];
    E2E1t = [currentEt Et(i)];
    E2E1b = [currentEb Eb(i)];
    E2E1p = [currentEp Ep(i)];
    E2E1c = [currentEc Ec(i)];
    E2E1m = [currentEm Em(i)];
    
    %CALCULO YIELD
    %Primero se calcula la sección eficaz media en el intervalo y
    %posterioermente el yield producido. Como se introduce el número de
    %protones por segundo en realidad se calcula el yield por unidad de
    %tiempo (s).
    
    % Water (full + simplified)
    Y_O16_C11(i) = rho_O16_A * trapz(E2E1, [O16_C11_F(E1)*1e-24/S1 O16_C11_F(E2)*1e-24/S2]);
    Y_O16_N13(i) = rho_O16_A * trapz(E2E1, [O16_N13_F(E1)*1e-24/S1 O16_N13_F(E2)*1e-24/S2]);
    Y_O16_O15(i) = rho_O16_A * trapz(E2E1, [O16_O15_F(E1)*1e-24/S1 O16_O15_F(E2)*1e-24/S2]);
    sigma_PG_C12_4w = 0.5 * (max(0,PG_C12_C12_4_F(E1)) + max(0,PG_C12_C12_4_F(E2)));
    sigma_PG_O16_4w = 0.5 * (max(0,PG_O16_C12_4_F(E1)) + max(0,PG_O16_C12_4_F(E2)));
    sigma_PG_N14_1w = 0.5 * (max(0,PG_N14_N14_1_F(E1)) + max(0,PG_N14_N14_1_F(E2)));
    sigma_PG_N14_2w = 0.5 * (max(0,PG_N14_N14_2_F(E1)) + max(0,PG_N14_N14_2_F(E2)));
    sigma_PG_O16_6w = 0.5 * (max(0,PG_O16_O16_6_F(E1)) + max(0,PG_O16_O16_6_F(E2)));
    sigma_C11_mean = 0.5 * (O16_C11_F(E1) + O16_C11_F(E2));
    sigma_N13_mean = 0.5 * (O16_N13_F(E1) + O16_N13_F(E2));
    sigma_O15_mean = 0.5 * (O16_O15_F(E1) + O16_O15_F(E2));
    Y_O16_C11s(i) = rho_O16_A * sigma_C11_mean * 1e-24 * dx;
    Y_O16_N13s(i) = rho_O16_A * sigma_N13_mean * 1e-24 * dx;
    Y_O16_O15s(i) = rho_O16_A * sigma_O15_mean * 1e-24 * dx;
    Y_PG_O16_C12_4w(i) = pps * rho_O16_A * sigma_PG_O16_4w * 1e-24 * dx;
    Y_PG_O16_O16_6w(i) = pps * rho_O16_A * sigma_PG_O16_6w * 1e-24 * dx;
    
    % Tissue (simplified only)
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
            
    
    Y_PG_trap(i) = pps * rho_O16_Ap * trapz(E2E1p, [PG_O16_O16_6_F(E1p)*1e-24/S1p PG_O16_O16_6_F(E2p)*1e-24/S2p]);
            
            

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
    
    % cartilage (simplified only)
    sigma_C11_meanc = 0.5 * (max(0,O16_C11_F(E1c)) + max(0,O16_C11_F(E2c)));
    sigma_N13_meanc = 0.5 * (max(0,O16_N13_F(E1c)) + max(0,O16_N13_F(E2c)));
    sigma_O15_meanc = 0.5 * (max(0,O16_O15_F(E1c)) + max(0,O16_O15_F(E2c)));
    sigma_N14_N13c = 0.5 * (max(0,N14_N13_F(E1c)) + max(0,N14_N13_F(E2c)));
    sigma_N14_O14c = 0.5 * (max(0,N14_O14_F(E1c)) + max(0,N14_O14_F(E2c)));
    sigma_N14_C11c = 0.5 * (max(0,N14_C11_F(E1c)) + max(0,N14_C11_F(E2c)));
    sigma_C12_C11c = 0.5 * (max(0,C12_C11_F(E1c)) + max(0,C12_C11_F(E2c)));
    sigma_C12_C10c = 0.5 * (max(0,C12_C10_F(E1c)) + max(0,C12_C10_F(E2c)));
    sigma_Ca44_Sc44c = 0.5 * (max(0,Ca44_Sc44_F(E1c)) + max(0,Ca44_Sc44_F(E2c)));
    sigma_PG_C12_4c = 0.5 * (max(0,PG_C12_C12_4_F(E1c)) + max(0,PG_C12_C12_4_F(E2c)));
    sigma_PG_O16_4c = 0.5 * (max(0,PG_O16_C12_4_F(E1c)) + max(0,PG_O16_C12_4_F(E2c)));
    sigma_PG_N14_1c = 0.5 * (max(0,PG_N14_N14_1_F(E1c)) + max(0,PG_N14_N14_1_F(E2c)));
    sigma_PG_N14_2c = 0.5 * (max(0,PG_N14_N14_2_F(E1c)) + max(0,PG_N14_N14_2_F(E2c)));
    sigma_PG_O16_6c = 0.5 * (max(0,PG_O16_O16_6_F(E1c)) + max(0,PG_O16_O16_6_F(E2c)));
    Y_O16_C11c(i) = pps *rho_O16_Ac * sigma_C11_meanc * 1e-24 * dx;
    Y_O16_N13c(i) = pps * rho_O16_Ac * sigma_N13_meanc * 1e-24 * dx;
    Y_O16_O15c(i) = pps * rho_O16_Ac * sigma_O15_meanc * 1e-24 * dx;
    Y_N14_N13c(i) = pps * rho_N14_Ac * sigma_N14_N13c * 1e-24 * dx;
    Y_N14_C11c(i) = pps * rho_N14_Ac * sigma_N14_C11c * 1e-24 * dx;! 
    Y_N14_O14c(i) = pps * rho_N14_Ac * sigma_N14_O14c * 1e-24 * dx;!
    Y_C12_C11c(i) = pps * rho_C12_Ac * sigma_C12_C11c * 1e-24 * dx;
    Y_C12_C10c(i) = pps * rho_C12_Ac * sigma_C12_C10c * 1e-24 * dx;
    Y_PG_C12_C12_4c(i) = pps * rho_C12_Ac * sigma_PG_C12_4c * 1e-24 * dx;
    Y_PG_O16_C12_4c(i) = pps * rho_O16_Ac * sigma_PG_O16_4c * 1e-24 * dx;
    Y_PG_N14_N14_1c(i) = pps * rho_N14_Ac * sigma_PG_N14_1c * 1e-24 * dx;
    Y_PG_N14_N14_2c(i) = pps * rho_N14_Ac * sigma_PG_N14_2c * 1e-24 * dx;
    Y_PG_O16_O16_6c(i) = pps * rho_O16_Ac * sigma_PG_O16_6c * 1e-24 * dx;
    
    
    %muscle
    sigma_C11_meanm = 0.5 * (max(0,O16_C11_F(E1m)) + max(0,O16_C11_F(E2m)));
    sigma_N13_meanm = 0.5 * (max(0,O16_N13_F(E1m)) + max(0,O16_N13_F(E2m)));
    sigma_O15_meanm = 0.5 * (max(0,O16_O15_F(E1m)) + max(0,O16_O15_F(E2m)));
    sigma_N14_N13m = 0.5 * (max(0,N14_N13_F(E1m)) + max(0,N14_N13_F(E2m)));
    sigma_N14_O14m = 0.5 * (max(0,N14_O14_F(E1m)) + max(0,N14_O14_F(E2m)));
    sigma_N14_C11m = 0.5 * (max(0,N14_C11_F(E1m)) + max(0,N14_C11_F(E2m)));
    sigma_C12_C11m = 0.5 * (max(0,C12_C11_F(E1m)) + max(0,C12_C11_F(E2m)));
    sigma_C12_C10m = 0.5 * (max(0,C12_C10_F(E1m)) + max(0,C12_C10_F(E2m)));
    sigma_Ca44_Sc44m = 0.5 * (max(0,Ca44_Sc44_F(E1m)) + max(0,Ca44_Sc44_F(E2m)));
    sigma_PG_C12_4m = 0.5 * (max(0,PG_C12_C12_4_F(E1m)) + max(0,PG_C12_C12_4_F(E2m)));
    sigma_PG_O16_4m = 0.5 * (max(0,PG_O16_C12_4_F(E1m)) + max(0,PG_O16_C12_4_F(E2m)));
    sigma_PG_N14_1m = 0.5 * (max(0,PG_N14_N14_1_F(E1m)) + max(0,PG_N14_N14_1_F(E2m)));
    sigma_PG_N14_2m = 0.5 * (max(0,PG_N14_N14_2_F(E1m)) + max(0,PG_N14_N14_2_F(E2m)));
    sigma_PG_O16_6m = 0.5 * (max(0,PG_O16_O16_6_F(E1m)) + max(0,PG_O16_O16_6_F(E2m)));
    Y_O16_C11m(i) = pps * rho_O16_Am * sigma_C11_meanm * 1e-24 * dx;
    Y_O16_N13m(i) = pps * rho_O16_Am * sigma_N13_meanm * 1e-24 * dx;
    Y_O16_O15m(i) = pps * rho_O16_Am * sigma_O15_meanm * 1e-24 * dx;
    Y_N14_N13m(i) = pps * rho_N14_Am * sigma_N14_N13m * 1e-24 * dx;
    Y_N14_C11m(i) = pps * rho_N14_Am * sigma_N14_C11m * 1e-24 * dx;! 
    Y_N14_O14m(i) = pps * rho_N14_Am * sigma_N14_O14m * 1e-24 * dx;!
    Y_C12_C11m(i) = pps * rho_C12_Am * sigma_C12_C11m * 1e-24 * dx;
    Y_C12_C10m(i) = pps * rho_C12_Am * sigma_C12_C10m * 1e-24 * dx;
    Y_PG_C12_C12_4m(i) = pps * rho_C12_Am * sigma_PG_C12_4m * 1e-24 * dx;
    Y_PG_O16_C12_4m(i) = pps * rho_O16_Am * sigma_PG_O16_4m * 1e-24 * dx;
    Y_PG_N14_N14_1m(i) = pps * rho_N14_Am * sigma_PG_N14_1m * 1e-24 * dx;
    Y_PG_N14_N14_2m(i) = pps * rho_N14_Am * sigma_PG_N14_2m * 1e-24 * dx;
    Y_PG_O16_O16_6m(i) = pps * rho_O16_Am * sigma_PG_O16_6m * 1e-24 * dx;
    
end

%% Create plots de emisores beta+

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
title('Yields of Water (per incoming proton)');
hold on
legend('C11','N13','O15','Location', 'Southeast');
ylabel('Yield');
plot(x,Y_O16_C11s,'k'); hold on
plot(x,Y_O16_N13s,'c')
plot(x,Y_O16_O15s,'m')
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
ylabel('Yield /mm');
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
ylabel('Yield/mm');
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
ylabel('Yield/mm');
set(gca,'FontSize',14)

% Figure in PMMA
figure
%subplot(2,1,1)
%plot(x,Ep)
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
plot(x,Y_C11p+Y_O15p+Y_N13p+Y_C10p,'b');
%plot(x,Y_C11p,'k');
%plot(x,Y_O15p,'m'); 
%plot(x,Y_N13p,'c');
%plot(x,Y_C10p,'y');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
%legend('Total','C11','O15','N13','C10','Location', 'northwest');
legend('Total','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 (max(Y_C11p+Y_O15p+Y_N13p+Y_C10p)+0.8*max(Y_C11p+Y_O15p+Y_N13p+Y_C10p))]);
xlabel('Depth (cm)');
ylabel('Yield /mm');
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
Y_C11p_norm = Y_C11p./(max(Y_C11p));
Y_CO_norm = (Y_C11p+Y_O15p)./max(Y_C11p+Y_O15p);
Y_tot_p=Y_C11p+Y_N13p+Y_O15p+Y_C10p;
plot(x,(Y_CO_norm),'b'); hold on
%plot(x,Y_CO_norm,'k');
% plot(x,Y_N13p,'c')
% plot(x,Y_O15p,'m')
% plot(x,Y_C10p,'y');
% plot(x,Y_tot_p,'.k');
% legend('C11','N13','O15','C10','Total','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 max(Y_CO_norm) + 1]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
set(gca,'FontSize',14)


% Figure in cartilage
figure
%subplot(2,1,1)
%plot(x,Ec)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Cartilage');
hold on
plot(x,100*Ddepc)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(100*Ddepc)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11c = Y_O16_C11c + Y_C12_C11c + Y_N14_C11c;
Y_N13c = Y_O16_N13c + Y_N14_N13c;
Y_O15c = Y_O16_O15c;
Y_C10c = Y_C12_C10c;
plot(x,Y_C11c,'k'); hold on
plot(x,Y_N13c,'c')
plot(x,Y_O15c,'m')
plot(x,Y_C10c,'y');
legend('C11','N13','O15','C10','Location', 'northwest');
[f,g]=min(Ddepc);
axis([0 (ceil(g*dx)) 0 max(max(Y_C11c,Y_O15c)+0.2*max(Y_C11c,Y_O15c))]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
set(gca,'FontSize',14)

% Figure in muscle
figure
%subplot(2,1,1)
%plot(x,Em)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Muscle');
hold on
plot(x,100*Ddepm)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(100*Ddepm)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11m = Y_O16_C11m + Y_C12_C11m + Y_N14_C11m;
Y_N13m = Y_O16_N13m + Y_N14_N13m;
Y_O15m = Y_O16_O15m;
Y_C10m = Y_C12_C10m;
plot(x,Y_C11m,'k'); hold on
plot(x,Y_N13m,'c')
plot(x,Y_O15m,'m')
plot(x,Y_C10m,'y');
legend('C11','N13','O15','C10','Location', 'northwest');
[f,g]=min(Ddepm);
axis([0 (ceil(g*dx)) 0 max(max(Y_C11m,Y_O15m)+0.2*max(Y_C11m,Y_O15m))]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
set(gca,'FontSize',14)

%% Figuras del Yield production de los PG

%PG Agua
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('WATER');
hold on
plot(x,100*Ddep)
legend('Dose')
set(gca,'FontSize',14)
axis([0 4 0 (max(100*Ddep)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4 = Y_PG_O16_C12_4w;
Y_6 = Y_PG_O16_O16_6w;
plot(x,Y_4,'r'); hold on
plot(x,Y_6,'m');
plot(x,Y_4+Y_6,'b');
% plot(x,Y_N13p,'c')
% plot(x,Y_O15p,'m')
% plot(x,Y_C10p,'y');
legend('4.44 MeV', '6.13 MeV','total','Location', 'northwest');
[f,g]=min(Ddep);
%axis([0 (ceil(g*dx)) 0 (max(Y_4+Y_6)+0.2*max(Y_4+Y_6))]);
axis([0 4 0 (max(Y_4+Y_6)+0.2*max(Y_4+Y_6))]);
xlabel('Depth (cm)');
ylabel('Isotopes/mm');
set(gca,'FontSize',14)

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
plot(x,Y_1p,'y');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('4.44 MeV', '6.13 MeV','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)) 0 (max(Y_4p)+0.2*max(Y_4p))]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
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
ylabel('Yield/mm');
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
ylabel('Yield/mm');
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
ylabel('Yield/mm');
set(gca,'FontSize',14)

%PG Cartilage
figure
%subplot(2,1,1)
%plot(x,Ec)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Cartilage');
hold on
plot(x,100*Ddepc)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepc)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4c = Y_PG_C12_C12_4c+Y_PG_O16_C12_4c;
Y_1c = Y_PG_N14_N14_1c;
Y_2c = Y_PG_N14_N14_2c;
Y_6c = Y_PG_O16_O16_6c;
plot(x,Y_6c,'b'); hold on
plot(x,Y_4c,'m');
plot(x,Y_2c,'c');
plot(x,Y_1c,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepc);
axis([0 (ceil(g*dx)) 0 (max(Y_4c)+0.2*max(Y_4c))]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
set(gca,'FontSize',14)

%PG Muscle
figure
%subplot(2,1,1)
%plot(x,Em)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Muscle');
hold on
plot(x,100*Ddepm)
legend('Dose')
set(gca,'FontSize',14)
axis([0 15 0 (max(100*Ddepm)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4m = Y_PG_C12_C12_4m+Y_PG_O16_C12_4m;
Y_1m = Y_PG_N14_N14_1m;
Y_2m = Y_PG_N14_N14_2m;
Y_6m = Y_PG_O16_O16_6m;
plot(x,Y_6m,'b'); hold on
plot(x,Y_4m,'m');
plot(x,Y_2m,'c');
plot(x,Y_1m,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepm);
axis([0 (ceil(g*dx)) 0 (max(Y_4m)+0.2*max(Y_4m))]);
xlabel('Depth (cm)');
ylabel('Yield/mm');
set(gca,'FontSize',14)

%% PET activity (t) para un cierto número de protones que llegan simultaneamente
%Valores para los que se pretende calcular la actividad
calcTimes = [0 60 1000 3600]; % s

%Definimos las variables
act_C11 = zeros(numel(calcTimes), numel(x)); %water
act_N13 = zeros(numel(calcTimes), numel(x));
act_O15 = zeros(numel(calcTimes), numel(x));

act_C11t = zeros(numel(calcTimes), numel(x)); %tissue
act_C10t = zeros(numel(calcTimes), numel(x));
act_N13t = zeros(numel(calcTimes), numel(x));
act_O15t= zeros(numel(calcTimes), numel(x));

act_C11a = zeros(numel(calcTimes), numel(x)); %adipose
act_C10a = zeros(numel(calcTimes), numel(x));
act_N13a = zeros(numel(calcTimes), numel(x));
act_O15a= zeros(numel(calcTimes), numel(x));

act_C11b = zeros(numel(calcTimes), numel(x)); %bone
act_C10b = zeros(numel(calcTimes), numel(x));
act_N13b = zeros(numel(calcTimes), numel(x));
act_O15b = zeros(numel(calcTimes), numel(x));
act_Sc44b = zeros(numel(calcTimes), numel(x));

act_C11p = zeros(numel(calcTimes), numel(x)); %PMMA
act_C10p = zeros(numel(calcTimes), numel(x));
act_N13p = zeros(numel(calcTimes), numel(x));
act_O15p = zeros(numel(calcTimes), numel(x));

act_C11c = zeros(numel(calcTimes), numel(x)); %cartilage
act_C10c = zeros(numel(calcTimes), numel(x));
act_N13c = zeros(numel(calcTimes), numel(x));
act_O15c= zeros(numel(calcTimes), numel(x));

act_C11m = zeros(numel(calcTimes), numel(x)); %muscle
act_C10m = zeros(numel(calcTimes), numel(x));
act_N13m = zeros(numel(calcTimes), numel(x));
act_O15m= zeros(numel(calcTimes), numel(x));

%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % Water
    act_C11(i,:) = deltat * landa_C11 .* Y_O16_C11s .* exp(- landa_C11 * calcTimes(i));
    act_N13(i,:) = deltat * landa_N13 .* Y_O16_N13s .* exp(- landa_N13 * calcTimes(i));
    act_O15(i,:) = deltat * landa_O15 .* Y_O16_O15s .* exp(- landa_O15 * calcTimes(i));
    
    % Tissue
    act_C11t(i,:) = deltat * landa_C11 .* Y_C11t .* exp(- landa_C11 * calcTimes(i));
    act_C10t(i,:) = deltat * landa_C10 .* Y_C10t .* exp(- landa_C10 * calcTimes(i));    
    act_N13t(i,:) = deltat * landa_N13 .* Y_N13t .* exp(- landa_N13 * calcTimes(i));
    act_O15t(i,:) = deltat * landa_O15 .* Y_O15t .* exp(- landa_O15 * calcTimes(i));  
    
     % Adipose
    act_C11a(i,:) = deltat * landa_C11 .* Y_C11a .* exp(- landa_C11 * calcTimes(i));
    act_C10a(i,:) = deltat * landa_C10 .* Y_C10a .* exp(- landa_C10 * calcTimes(i));    
    act_N13a(i,:) = deltat * landa_N13 .* Y_N13a .* exp(- landa_N13 * calcTimes(i));
    act_O15a(i,:) = deltat * landa_O15 .* Y_O15a .* exp(- landa_O15 * calcTimes(i)); 
    
     % Bone
    act_C11b(i,:) = deltat * landa_C11 .* Y_C11b .* exp(- landa_C11 * calcTimes(i));
    act_C10b(i,:) = deltat * landa_C10 .* Y_C10b .* exp(- landa_C10 * calcTimes(i));    
    act_N13b(i,:) = deltat * landa_N13 .* Y_N13b .* exp(- landa_N13 * calcTimes(i));
    act_O15b(i,:) = deltat * landa_O15 .* Y_O15b .* exp(- landa_O15 * calcTimes(i)); 
    act_Sc44b(i,:) = deltat * landa_Sc44 .* Y_Sc44b .* exp(- landa_Sc44 * calcTimes(i));
    
     % PMMA
    act_C11p(i,:) = deltat * landa_C11 .* Y_C11p .* exp(- landa_C11 * calcTimes(i));
    act_C10p(i,:) = deltat * landa_C10 .* Y_C10p .* exp(- landa_C10 * calcTimes(i));    
    act_N13p(i,:) = deltat * landa_N13 .* Y_N13p .* exp(- landa_N13 * calcTimes(i));
    act_O15p(i,:) = deltat * landa_O15 .* Y_O15p .* exp(- landa_O15 * calcTimes(i)); 
    
      % Cartilage
    act_C11c(i,:) = deltat * landa_C11 .* Y_C11c .* exp(- landa_C11 * calcTimes(i));
    act_C10c(i,:) = deltat * landa_C10 .* Y_C10c .* exp(- landa_C10 * calcTimes(i));    
    act_N13c(i,:) = deltat * landa_N13 .* Y_N13c .* exp(- landa_N13 * calcTimes(i));
    act_O15c(i,:) = deltat * landa_O15 .* Y_O15c .* exp(- landa_O15 * calcTimes(i));
    
      % Muscle
    act_C11m(i,:) = deltat * landa_C11 .* Y_C11m .* exp(- landa_C11 * calcTimes(i));
    act_C10m(i,:) = deltat * landa_C10 .* Y_C10m .* exp(- landa_C10 * calcTimes(i));    
    act_N13m(i,:) = deltat * landa_N13 .* Y_N13m .* exp(- landa_N13 * calcTimes(i));
    act_O15m(i,:) = deltat * landa_O15 .* Y_O15m .* exp(- landa_O15 * calcTimes(i)); 
    
end
act_total = (act_C11 + act_N13 + act_O15);
act_totalt = (act_C11t + act_C10t + act_N13t + act_O15t);
act_totala = (act_C11a + act_C10a + act_N13a + act_O15a);
act_totalb = (act_C11b + act_C10b + act_N13b + act_O15b+act_Sc44b);
act_totalp = (act_C11p + act_C10p + act_N13p + act_O15p);
act_totalc = (act_C11c + act_C10c + act_N13c + act_O15c);
act_totalm = (act_C11m + act_C10m + act_N13m + act_O15m);

%Representación gráfica de los 4 tiempos calculados
F=figure;
G=figure;
H=figure;
I=figure;
J=figure;
K=figure;
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
    
    figure(J)
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totalc(i,:),'b')
    hold on
    plot(x, act_C11c(i,:),'r')
    plot(x, act_N13c(i,:),'y')
    plot(x, act_O15c(i,:),'c')
    plot(x, act_C10c(i,:),'g')
    title(sprintf('Activity at t=%i s',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
    axis([0 5 0 max(act_totalc(i,:))+0.2*max(act_totalc(i,:))]);   
    yyaxis right
    grid on
    plot(x,100*Ddepc);
    ylabel('-dE/dx')
    [f,g]=min(Ddepc);
    axis([0 (ceil(g*dx)) 0 max(100*Ddepc)+0.2*max(100*Ddepc)]);
    
    figure(K)
    subplot(2,2,i)
    yyaxis left
    plot(x, act_totalm(i,:),'b')
    hold on
    plot(x, act_C11m(i,:),'r')
    plot(x, act_N13m(i,:),'y')
    plot(x, act_O15m(i,:),'c')
    plot(x, act_C10m(i,:),'g')
    title(sprintf('Activity at t=%i s',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
    axis([0 5 0 max(act_totalm(i,:))+0.2*max(act_totalm(i,:))]);   
    yyaxis right
    grid on
    plot(x,100*Ddepm);
    ylabel('-dE/dx')
    [f,g]=min(Ddepm);
    axis([0 (ceil(g*dx)) 0 max(100*Ddepm)+0.2*max(100*Ddepm)]);
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
%Aquí se calcula la actividad en función del tiempo cuando el tiempo de
%irradiación es superior al delta de t, es decir, un caso clínico real
%usando un ciclotrón (no está preparado todavía para la simualción de
%pulsos de un sincrotrón). Los parámetros de tiempo de irradiación y de
%decaimiento se introducen al principio

c=1;
temp_C10t=zeros(t+1,numel(x));temp_C10a=zeros(t+1,numel(x));
temp_C11t=zeros(t+1,numel(x));temp_C11a=zeros(t+1,numel(x));
temp_N13t=zeros(t+1,numel(x));temp_N13a=zeros(t+1,numel(x));
temp_O15t=zeros(t+1,numel(x));temp_O15a=zeros(t+1,numel(x));
temp_C10p=zeros(t+1,numel(x));temp_C10b=zeros(t+1,numel(x));
temp_C11p=zeros(t+1,numel(x));temp_C11b=zeros(t+1,numel(x));
temp_N13p=zeros(t+1,numel(x));temp_N13b=zeros(t+1,numel(x));
temp_O15p=zeros(t+1,numel(x));temp_O15b=zeros(t+1,numel(x));
temp_C10c=zeros(t+1,numel(x));temp_C10m=zeros(t+1,numel(x));
temp_C11c=zeros(t+1,numel(x));temp_C11m=zeros(t+1,numel(x));
temp_N13c=zeros(t+1,numel(x));temp_N13m=zeros(t+1,numel(x));
temp_O15c=zeros(t+1,numel(x));temp_O15m=zeros(t+1,numel(x));

temp_totalp2=zeros(t+1);
temp_totalt2=zeros(t+1);
temp_totalb2=zeros(t+1);
temp_totala2=zeros(t+1);
temp_totalc2=zeros(t+1);
temp_totalm2=zeros(t+1);
temp_parcC11t=zeros(t+1);temp_parcC11p=zeros(t+1);temp_parcC11a=zeros(t+1);temp_parcC11b=zeros(t+1);temp_parcC11c=zeros(t+1);temp_parcC11m=zeros(t+1);
temp_parcO15t=zeros(t+1);temp_parcO15p=zeros(t+1);temp_parcO15a=zeros(t+1);temp_parcO15b=zeros(t+1);temp_parcO15c=zeros(t+1);temp_parcO15m=zeros(t+1);
temp_parcN13t=zeros(t+1);temp_parcN13p=zeros(t+1);temp_parcN13a=zeros(t+1);temp_parcN13b=zeros(t+1);temp_parcN13c=zeros(t+1);temp_parcN13m=zeros(t+1);
temp_parcC10t=zeros(t+1);temp_parcC10p=zeros(t+1);temp_parcC10a=zeros(t+1);temp_parcC10b=zeros(t+1);temp_parcC10c=zeros(t+1);temp_parcC10m=zeros(t+1);

%CALCULO
%Para hacer el calculo suponemos que cada intervalod de tiempo llegan un
%número de protones que se frenan inmediatamente y genera los isótopos
%correspondientes. En cada iteración se calcula el número de los isótopos
%de las iteraciones anteriores que se ha desintegrado y además si estamos
%en el tiempo de irradiación se suma la contribución de los protones que
%llegan.
for i=1:t+1
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
            
            temp_C10c(i,:)=temp_C10c(i,:)+deltat * Y_C10c.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11c(i,:)=temp_C11c(i,:)+deltat * Y_C11c.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13c(i,:)=temp_N13c(i,:)+deltat * Y_N13c.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15c(i,:)=temp_O15c(i,:)+deltat * Y_O15c.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
            
            temp_C10m(i,:)=temp_C10m(i,:)+deltat * Y_C10m.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1));
            temp_C11m(i,:)=temp_C11m(i,:)+deltat * Y_C11m.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1));
            temp_N13m(i,:)=temp_N13m(i,:)+deltat * Y_N13m.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));
            temp_O15m(i,:)=temp_O15m(i,:)+deltat * Y_O15m.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1));
            
        end
    else
            for j=1:a
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
            
            temp_C10c(i,:)=temp_C10c(i,:)+deltat *Y_C10c.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11c(i,:)=temp_C11c(i,:)+deltat *Y_C11c.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13c(i,:)=temp_N13c(i,:)+deltat *Y_N13c.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15c(i,:)=temp_O15c(i,:)+deltat *Y_O15c.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c)); 
            
            temp_C10m(i,:)=temp_C10m(i,:)+deltat *Y_C10m.*landa_C10.*exp(-landa_C10*(deltat*j-deltat*1+deltat*c));
            temp_C11m(i,:)=temp_C11m(i,:)+deltat *Y_C11m.*landa_C11.*exp(-landa_C11*(deltat*j-deltat*1+deltat*c));
            temp_N13m(i,:)=temp_N13m(i,:)+deltat *Y_N13m.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));
            temp_O15m(i,:)=temp_O15m(i,:)+deltat *Y_O15m.*landa_O15.*exp(-landa_O15*(deltat*j-deltat*1+deltat*c)); 
            
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

    temp_totalc2(i)=sum(temp_C10c(i,:))+sum(temp_C11c(i,:))+sum(temp_N13c(i,:))+sum(temp_O15c(i,:));
    temp_parcC11c(i)=sum(temp_C11c(i,:));
    temp_parcN13c(i)=sum(temp_N13c(i,:));
    temp_parcC10c(i)=sum(temp_C10c(i,:));
    temp_parcO15c(i)=sum(temp_O15c(i,:));
    
    temp_totalm2(i)=sum(temp_C10m(i,:))+sum(temp_C11m(i,:))+sum(temp_N13m(i,:))+sum(temp_O15m(i,:));
    temp_parcC11m(i)=sum(temp_C11m(i,:));
    temp_parcN13m(i)=sum(temp_N13m(i,:));
    temp_parcC10m(i)=sum(temp_C10m(i,:));
    temp_parcO15m(i)=sum(temp_O15m(i,:));
    
end
    temp_totalt=temp_C10t+temp_C11t+temp_N13t+temp_O15t;
    temp_totalp=temp_C10p+temp_C11t+temp_N13p+temp_O15p;
    temp_totala=temp_C10a+temp_C11a+temp_N13a+temp_O15a;
    temp_totalb=temp_C10b+temp_C11b+temp_N13b+temp_O15b;
    temp_totalc=temp_C10c+temp_C11c+temp_N13c+temp_O15c;
    temp_totalm=temp_C10m+temp_C11m+temp_N13m+temp_O15m;



    %% Dibuja la Actividad total en función del tiempo.
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totalt2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15t(:,1));
    plot(T*deltat,temp_parcC11t(:,1));
    plot(T*deltat,temp_parcN13t(:,1));
    plot(T*deltat,temp_parcC10t(:,1));
    title('Actividad total en a lo largo del tiempo PIEL');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters/s /proton/mm');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalt2(:,1))+0.2*max(temp_totalt2(:,1)))]); 
    
    figure;
    plot(T*deltat,(temp_totalp2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15p(:,1));
    plot(T*deltat,temp_parcC11p(:,1));
    plot(T*deltat,temp_parcN13p(:,1));
    plot(T*deltat,temp_parcC10p(:,1));
    title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters /proton/mm ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(temp_totalp2(:,1))+0.2*max(temp_totalp2(:,1)))]);
    
    figure;
    plot(T*deltat,(temp_totala2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15a(:,1));
    plot(T*deltat,temp_parcC11a(:,1));
    plot(T*deltat,temp_parcN13a(:,1));
    plot(T*deltat,temp_parcC10a(:,1));
    title('Actividad total en a lo largo del tiempo ADIPOSE');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters /proton/mm  ');
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
    ylabel('Beta+ emitters /proton/mm ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(temp_totalb2(:,1))+0.2*max(temp_totalb2(:,1)))]);
    
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totalc2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15c(:,1));
    plot(T*deltat,temp_parcC11c(:,1));
    plot(T*deltat,temp_parcN13c(:,1));
    plot(T*deltat,temp_parcC10c(:,1));
    title('Actividad total en a lo largo del tiempo CARTÍLAGO');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters /proton/mm  ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalc2(:,1))+0.2*max(temp_totalc2(:,1)))]); 
    
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totalm2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15m(:,1));
    plot(T*deltat,temp_parcC11m(:,1));
    plot(T*deltat,temp_parcN13m(:,1));
    plot(T*deltat,temp_parcC10m(:,1));
    title('Actividad total en a lo largo del tiempo MÚSCULO');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters /proton/mm  ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalm2(:,1))+0.2*max(temp_totalm2(:,1)))]); 
    
%% Actividad total en función de tiempo
%Igual que el apartado anterior pero con otra fórmula. Es una comprobación
%de que estaba bien hecho.
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
    plot(T*deltat,(int_totalp2(:,1)));
    hold on
    plot(T*deltat,int_parcO15p(:,1))
    plot(T*deltat,int_parcC11p(:,1))
    plot(T*deltat,int_parcC10p(:,1))
    plot(T*deltat,int_parcN13p(:,1))
    title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)')
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','C10','N13','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    axis([0 t*deltat 0 (max(int_totalp2(:,1))+0.2*max(int_totalp2(:,1)))]);

    
    
%% Actividad durante Bean-On
%Aquí calculamos la cantidad de pares de 511 keV que se producen durante la
%irradiación y los que ocurren después en función del espacio, hasta el valor del tiempo en
%segundos que hay de la variable de abajo
tt=240;

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

for i=1:tt
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


%% Actividad/Grey 
    %en el pico de Bragg 1cm3
Q = figure;
L = figure;
M = figure;
vol = 0.001; %cm3 
Ec_J = 1.6e-19*10^6*max(Ddepc);
rho_cartilage_kgcm3 = rho_cartilage*10^-3;
Ddepc_Gy_max=(Ec_J/(rho_cartilage_kgcm3 * vol));
Ddepc_Gy = Ddepc.*(Ddepc_Gy_max/max(Ddepc))*pps;
Em_J = 1.6e-19*10^6.*max(Ddepm);
rho_muscle_kgcm3 = rho_muscle*10^-3;
Ddepm_Gy_max=(Em_J/(rho_muscle_kgcm3 * vol));
Ddepm_Gy = Ddepm.*(Ddepm_Gy_max/max(Ddepm))*pps;

for i=1:numel(calcTimes)
    % cartilage
    A_cartilage = act_totalc(i,:) ./(Ec_J/(rho_cartilage_kgcm3 * vol));
    figure(Q)
    subplot(2,2,i)
    yyaxis right
    plot(x,Ddepc_Gy,'r')
    ylabel('Dose (Gy)')
    axis([0 3.5 0 (max(Ddepc_Gy)+max(Ddepc_Gy/10))])
    hold on
    yyaxis left
    plot(x,A_cartilage,'b')
    title(sprintf('Activity at t=%i s (cartilage) ',calcTimes(i)));
    xlabel('depth (cm)')
    ylabel('activity / dose V (Bq / Gy)')
    
    
    % muscle
    A_muscle = act_totalm(i,:) ./(Em_J/(rho_muscle_kgcm3 * vol));
    figure(L)
    subplot(2,2,i)
    yyaxis right
    plot(x,Ddepm_Gy,'r')
    ylabel('Dose (Gy)')
    axis([0 3.5 0 (max(Ddepm_Gy)+max(Ddepm_Gy/10))])
    hold on
    yyaxis left
    plot(x,A_muscle,'b')
    title(sprintf('Activity at t=%i s (muscle)',calcTimes(i)));
    xlabel('depth (cm)')
    ylabel('activity / dose (Bq / Gy)')

    % both
    figure(M)
    subplot(2,2,i)
    yyaxis right
    plot(x,Ddepm_Gy,'-r')
    hold on
    plot(x,Ddepc_Gy,'--r')
    axis([0 3.5 0 (max(Ddepc_Gy)+max(Ddepc_Gy/10))])
    ylabel('Dose (Gy)')
    hold on
    yyaxis left
    plot(x,A_muscle,'b')
    hold on
    plot(x,A_cartilage,'--b')
    title(sprintf('Activity at t=%i s (c-m)',calcTimes(i)));
    legend('muscle dose','cartilage dose','muscle','cartilage','Location','southwest')
    xlabel('depth (cm)')
    ylabel('activity / dose (Bq / Gy)')
    
    

end

%% Yield/Gy


% Figure in cartilage
figure
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (Gy)')
title('Cartilage');
hold on
Ddepc_Gy = Ddepc.*(Ddepc_Gy_max/max(Ddepc));
plot(x,Ddepc_Gy)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(Ddepc_Gy)+max(Ddepc_Gy/10))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11c_Gy = Y_C11c./(Ec_J/(rho_cartilage_kgcm3 * vol));
Y_N13c_Gy = Y_N13c./(Ec_J/(rho_cartilage_kgcm3 * vol));
Y_O15c_Gy = Y_O15c./(Ec_J/(rho_cartilage_kgcm3 * vol));
Y_C10c_Gy = Y_C10c./(Ec_J/(rho_cartilage_kgcm3 * vol));
plot(x,Y_C11c_Gy,'k'); hold on
plot(x,Y_N13c_Gy,'c')
plot(x,Y_O15c_Gy,'m')
plot(x,Y_C10c_Gy,'y');
legend('C11','N13','O15','C10','Location', 'northeast');
ylabel('\beta^+ emitters/Gy')
xlabel('Depth')
title('cartilage')

% Figure in cartilage
figure
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (Gy)')
title('Muscle');
hold on
Ddepm_Gy = Ddepc.*(Ddepc_Gy_max/max(Ddepc));
plot(x,Ddepm_Gy)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(Ddepm_Gy)+max(Ddepm_Gy/10))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11m_Gy = Y_C11m./(Em_J/(rho_muscle_kgcm3 * vol));
Y_N13m_Gy = Y_N13m./(Em_J/(rho_muscle_kgcm3 * vol));
Y_O15m_Gy = Y_O15m./(Em_J/(rho_muscle_kgcm3 * vol));
Y_C10m_Gy = Y_C10m./(Em_J/(rho_muscle_kgcm3 * vol));
plot(x,Y_C11m_Gy,'k'); hold on
plot(x,Y_N13m_Gy,'c')
plot(x,Y_O15m_Gy,'m')
plot(x,Y_C10m_Gy,'y');
legend('C11','N13','O15','C10','Location', 'northeast');
ylabel('\beta^+ emitters/Gy')
xlabel('Depth')
title('muscle')


