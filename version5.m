% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es

%% Par?metros a modificar:
clear all;close all;

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control1.mat');

%PARAMETROS
dx=0.02;      %Paso del intervalo (cm)
xref=20;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=140;        %Energía inicial del haz
deltat=1;      %Inervalo de tiempo de las simulaciones
a=120/deltat;  %Tiempo de irradación del haz (s)
t=900/deltat;  %Tiempo total de la simulación
tt=240/deltat; %Tiempo de recogida de datos total
pps=1; %protones/segundo
MeVJ=1.6e-13;
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

%Densidades Atomicas (A falta de las masas molares usamos
%sum(Comp_bone.*W_ele)) CORREGIR
!rho_tissue_A = 4.6243E+22; % atoms/cm3
rho_tissue_A = AvNmbr*rho_tissue/sum(Comp_tissue.*W_ele);  % atoms/cm3
rho_bone_A = AvNmbr*rho_bone/sum(Comp_bone.*W_ele);  % atoms/cm3
rho_adipose_A = AvNmbr*rho_adipose/sum(Comp_adipose.*W_ele);  % atoms/cm3
rho_PMMA_A = AvNmbr*rho_PMMA/sum(Comp_PMMA.*W_ele);  % atoms/cm3
%rho_PMMA_A = AvNmbr*rho_PMMA/PMMA_Molar;  % atoms/cm3
rho_w_A =  rho_w * AvNmbr / waterMolecularWeight; % molecules / cm3
ZnAtomicWeight = 65.38; % g/mol
rho_Zn_A = Zn_fraction * rho_Zn * AvNmbr / ZnAtomicWeight; % molecules / cm3

%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
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

%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
E = nan(size(x));
Et = nan(size(x));
Eb = nan(size(x));
Ea = nan(size(x));
Ep = nan(size(x));

%Energia depositada por material
Ddep = zeros(size(E));
Ddept = zeros(size(E));
Ddepb = zeros(size(E));
Ddepa = zeros(size(E));
Ddepp = zeros(size(E));

%Energia actual de cada material
currentE = E0;
currentEt = E0;
currentEb = E0;
currentEa = E0;
currentEp = E0;

%CREACION DE VECTORES DE YIELD
% In water (full + simplified versions)
Y_O16_C11 = zeros(size(x));
Y_O16_N13 = zeros(size(x));
Y_O16_O15 = zeros(size(x));
Y_O16_C11s = zeros(size(x));
Y_O16_N13s = zeros(size(x));
Y_O16_O15s = zeros(size(x));
Y_O18_F18w = zeros(size(x));
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
    
    %Multiplicamospor la densidad al poder de frenado
    S1 = (S_w*rho_w); % MeV/cm
    S1t = (S_t*rho_tissue); %*(1-Zn_fraction) + S_Znt*rho_Zn*Zn_fraction); % MeV/cm
    S1b = (S_b*rho_bone); % MeV/cm
    S1a = (S_a*rho_adipose); % MeV/cm
    S1p = (S_p*rho_PMMA); % MeV/cm
    
    %Calculamos la energía que pierde el protón al atravesar la lámina (se
    %supone que toda la pérdida se hace al final).
    deltaE = dx*S1; % MeV
    deltaEt = dx*S1t; % MeV
    deltaEb = dx*S1b; % MeV
    deltaEa = dx*S1a; %MeV
    deltaEp = dx*S1p; %MeV
    
    %Guardamos las energías para poder representar todo luego en función de
    %ella.
    E(i) = currentE; % MeV
    Et(i) = currentEt; % MeV
    Eb(i) = currentEb;  % MeV
    Ea(i) = currentEa; % MeV
    Ep(i) = currentEp; % MeV
    
    %Calculamos la energía con la que sale de la lámina que se usará en la
    %siguiente iteración.
    currentE = currentE - deltaE; % MeV
    currentEt = currentEt - deltaEt; % MeV
    currentEb = currentEb - deltaEb; % MeV
    currentEa = currentEa - deltaEa; %MeV
    currentEp = currentEp - deltaEp; %MeV
    
    %NO SE USA
    S_w2 = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Zn2 = max(0,1000*S_Zn_F(currentE*1000)); % MeV/(g/cm2)
    S_t2 = max(0,1000*S_t_F(currentEt*1000)); % MeV/(g/cm2)
    S_Zn2t = max(0,1000*S_Zn_F(currentEt*1000)); % MeV/(g/cm2)
    S_b2 = max(0,1000*S_b_F(currentEb*1000));% MeV/(g/cm2)
    S_a2 = max(0,1000*S_a_F(currentEa*1000));% MeV/(g/cm2)
    S_p2 = max(0,1000*S_p_F(currentEp*1000));% MeV/(g/cm2)
    S2 = (S_w2*rho_w); % MeV/cm
    S2t = (S_t2*rho_tissue); %*(1-Zn_fraction) + S_Zn2t*rho_Zn*Zn_fraction); % MeV/cm
    S2b = (S_b2*rho_bone); % MeV/cm3
    S2a = (S_a2*rho_adipose); % MeV/cm3
    S2p = (S_p2*rho_PMMA); % MeV/cm3
    
    %Guardamos la energía depositada en cada uno de los materiales.
    Ddep(i) = deltaE; % MeV
    Ddept(i) = deltaEt; % MeV
    Ddepb(i) = deltaEb;  % MeV
    Ddepa(i) = deltaEa; %MeV
    Ddepp(i) = deltaEp; %MeV
    
    %Creamos unas variables con los datos anteriores para introducir a
    %continuación en los Yield.
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
    %Estas de abajo solo se usan por el método trapezoidal que no es el que
    %estamos teniendo en cuenta.
    E2E1 = [currentE E(i)];
    E2E1t = [currentEt Et(i)];
    E2E1b = [currentEb Eb(i)];
    E2E1p = [currentEp Ep(i)];
    
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
    sigma_F18_mean = 0.5 * (O18_F18_F(E1) + O18_F18_F(E2));
    Y_O16_C11s(i) = rho_O16_A * sigma_C11_mean * 1e-24 * dx;
    Y_O16_N13s(i) = rho_O16_A * sigma_N13_mean * 1e-24 * dx;
    Y_O16_O15s(i) = rho_O16_A * sigma_O15_mean * 1e-24 * dx;
    Y_O18_F18w(i) = rho_O16_A * sigma_F18_mean * 1e-24 * dx;
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
    
    
end

%% Calculo de unidades Dosis

%Water
Edep=0;
for i=1:length(x)
    if Ddep(i)>0.2*max(Ddep)
        Edep=Edep+Ddep(i);
    end
end
Da_w=Edep*1.6e-13/rho_w;
Np_w=1/Da_w;
Con_w=Np_w*MeVJ;

%Tissue
Edept=0;
for i=1:length(x)
    if Ddept(i)>0.2*max(Ddept)
        Edept=Edept+Ddept(i);
    end
end
Da_t=Edept*1.6e-13/rho_tissue;
Np_t=1/Da_t;
Con_t=Np_t*MeVJ;

%Adipose
Edepa=0;
for i=1:length(x)
    if Ddepa(i)>0.2*max(Ddepa)
        Edepa=Edepa+Ddepa(i);
    end
end
Da_a=Edepa*1.6e-13/rho_adipose;
Np_a=1/Da_a;
Con_a=Np_a*MeVJ;

%Bone
Edepb=0;
for i=1:length(x)
    if Ddepb(i)>0.2*max(Ddepb)
        Edepb=Edepb+Ddepb(i);
    end
end
Da_b=Edepb*1.6e-13/rho_bone;
Np_b=1/Da_b;
Con_b=Np_b*MeVJ;

%PMMA
Edepp=0;
for i=1:length(x)
    if Ddepp(i)>0.2*max(Ddepp)
        Edepp=Edepp+Ddepp(i);
    end
end
Da_p=Edep*1.6e-13/rho_PMMA;
Np_p=1/Da_p;
Con_p=Np_p*MeVJ;

%% Create plots de emisores beta+

% Figure in water
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,E)
yyaxis right
xlabel('Depth (cm)');
%hold on
plot(x,Con_w*Ddep,'linewidth',2)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_w*Ddep)+0.2*max(Con_w*Ddep))]);
%subplot(2,1,2)
yyaxis left
title(' Water ');
hold on
ylabel('\beta^+ isotopes/Gy/mm');
Y_wt = Y_O16_C11s+Y_O16_N13s+Y_O16_O15s;
plot(x,Np_w/pps*Y_O16_C11s,'r'); hold on
plot(x,Np_w/pps*Y_O16_N13s,'c')
plot(x,Np_w/pps*Y_O16_O15s,'m')
plot(x,Np_w/pps*Y_wt,'k','linewidth',2);
legend('C11','N13','O15','Total','Location', 'northeastoutside');
set(gca,'FontSize',14)
grid on
[f,g]=min(Ddep);
axis([g*dx-1 g*dx+0.00*g*dx  0 max(Np_w/pps*Y_wt)+0.2*max(Np_w/pps*Y_wt)]);



% Figure in tisssue
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,Et)
yyaxis right
xlabel('Depth (cm)');
ylabel(' Dose (Gy/cm³)')
title('Tissue');
hold on
plot(x,Con_t*Ddept,'linewidth',2)
legend('Dose')
set(gca,'FontSize',14)
axis([0 17 0 (max(Con_t*Ddept)+0.2*max(Con_t*Ddept))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11t = Y_O16_C11t + Y_C12_C11t + Y_N14_C11t;
Y_N13t = Y_O16_N13t + Y_N14_N13t;
Y_O15t = Y_O16_O15t;
Y_C10t = Y_C12_C10t;
Y_tt = Y_C11t+Y_N13t+Y_O15t+Y_C10t;
plot(x,Np_t/pps*Y_C11t,'r'); hold on
plot(x,Np_t/pps*Y_N13t,'c')
plot(x,Np_t/pps*Y_O15t,'m')
plot(x,Np_t/pps*Y_C10t,'b');
plot(x,Np_t/pps*Y_tt,'k','linewidth',2);
legend('C11','N13','O15','C10','Total','Location', 'northeastoutside');
[f,g]=min(Ddept);
axis([g*dx-1 g*dx+0.00*g*dx 0 max(Np_t/pps*Y_tt)+0.2*max(Np_t/pps*Y_tt)]);
xlabel('Depth (cm)');
ylabel('\beta^+ isotopes/Gy/mm');
set(gca,'FontSize',14)

% Figure in adipose
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,Ea)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (Gy/cm³)')
title('Adipose');
hold on
plot(x,Con_a*Ddepa,'linewidth',2)
legend('Dose')
set(gca,'FontSize',14)
axis([0 20 0 (max(Con_a*Ddepa)+0.2*max(Con_a*Ddepa))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_C11a = Y_O16_C11a + Y_C12_C11a + Y_N14_C11a;
Y_N13a = Y_O16_N13a + Y_N14_N13a;
Y_O15a = Y_O16_O15a;
Y_C10a = Y_C12_C10a;
Y_ta = Y_C11a+Y_N13a+Y_O15a+Y_C10a;
plot(x,Np_a/pps*Y_C11a,'k'); hold on
plot(x,Np_a/pps*Y_N13a,'c')
plot(x,Np_a/pps*Y_O15a,'m')
plot(x,Np_a/pps*Y_C10a,'y');
plot(x,Np_a/pps*Y_ta,'k','linewidth',2);
legend('C11','N13','O15','C10','Location', 'northeastoutside');
[f,g]=min(Ddepa);
axis([g*dx-1 g*dx+0.00*g*dx 0 max(Np_a/pps*Y_ta)+0.2*max(Np_a/pps*Y_ta)]);
xlabel('Depth (cm)');
ylabel('\beta^+ isotopes/Gy/mm');
set(gca,'FontSize',14)

% Figure in bone
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (Gy/cm³)')
title('Bone');
hold on
plot(x,Con_b*Ddepb,'linewidth',2)
legend('Dose')
set(gca,'FontSize',14)
axis([0 11 0 (max(Con_b*Ddepb)+0.2*max(Con_b*Ddepb))]);
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
Y_tb = Y_C11b+Y_N13b+Y_O15b+Y_C10b+Y_Sc44b;
plot(x,Np_b/pps*Y_C11b,'k'); hold on
plot(x,Np_b/pps*Y_N13b,'c')
plot(x,Np_b/pps*Y_O15b,'m')
plot(x,Np_b/pps*Y_C10b,'y');
plot(x,Np_b/pps*Y_Sc44b,'g');
plot(x,Np_b/pps*Y_tb,'k','linewidth',2);
legend('C11','N13','O15','C10','Sc44','Total','Location', 'northeastoutside');
[f,g]=min(Ddepb);
axis([g*dx-1 g*dx+0.00*g*dx 0 max(Np_b/pps*Y_tb)+0.2*max(Np_b/pps*Y_tb)]);
xlabel('Depth (cm)');
ylabel('\beta^+ isotopes/Gy/mm');
set(gca,'FontSize',14)

% % Figure in PMMA
% figure
% %subplot(2,1,1) plot(x,Eb)
% yyaxis right
% xlabel('Depth (cm)');
% ylabel('Dose (Gy/cm³)')
% title('PMMA');
% hold on
% plot(x,Con_p*Ddepp)
% legend('Dose')
% set(gca,'FontSize',14)
% axis([0 15 0 (max(Con_p*Ddepp)+0.2*max(Con_p*Ddepp))]);
% %subplot(2,1,2)
% yyaxis left
% grid on
% %title('Yields of different species (per incoming proton)');
% hold on
% Y_C11p = Y_O16_C11p + Y_C12_C11p + Y_N14_C11p;
% Y_N13p = Y_O16_N13p + Y_N14_N13p;
% Y_O15p = Y_O16_O15p;
% Y_C10p = Y_C12_C10p;
% hold on
% %plot(x,Np_p/pps*Y_C11p,'b'); plot(x,Np_p/pps*Y_O15p,'m');
% %plot(x,Np_p/pps*Y_N13p,'c'); plot(x,Np_p/pps*Y_C10p,'y');
% %plot(x,Np_p/pps*Y_C11p+Np_p/pps*Y_O15p+Np_p/pps*Y_N13p+Np_p/pps*Y_C10p,'k');
% %plot(x,Y_N13p,'c')
% plot(x,Np_p/pps*Y_C11p+Np_p/pps*Y_O15p,'k')
% %plot(x,Y_O15p,'m') plot(x,Y_C10p,'y');
% %legend('C11','O15','N13','C10','Total','Location', 'northwest');
% %legend('Total','Location', 'northwest');
% legend('C11 + O15','Location', 'northwest');
% [f,g]=min(Ddepp);
% axis([0 (ceil(g*dx)) 0 (max(Np_p/pps*Y_C11p)+0.8*max(Np_p/pps*Y_C11p))]);
% xlabel('Depth (cm)');
% ylabel('C11+O15 \beta^+ isotopes/Gy/mm');
% set(gca,'FontSize',14)


%% Figuras del Yield production de los PG

%PG Agua
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (Gy/cm³)')
title('WATER');
hold on
plot(x,Con_w*Ddep)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_w*Ddep)+0.2*max(Con_w*Ddep))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4 = Y_PG_O16_C12_4w;
Y_6 = Y_PG_O16_O16_6w;
plot(x,Np_w/pps*Y_4,'r'); hold on
plot(x,Np_w/pps*Y_6,'m');
%plot(x,Np_w/pps*Y_4+Np_w/pps*Y_6,'b');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('4.44 MeV', '6.13 MeV','Location', 'northwest');
[f,g]=min(Ddep);
%axis([0 (ceil(g*dx)) 0 (max(Y_4+Y_6)+0.2*max(Y_4+Y_6))]);
axis([0 (ceil(g*dx)) 0 (max(Np_w/pps*Y_4)+0.2*max(Np_w/pps*Y_4))]);
xlabel('Depth (cm)');
ylabel('PG/Gy/mm');
set(gca,'FontSize',14)

%PG PMMA
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
title('PMMA');
hold on
plot(x,Con_p*Ddepp)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_p*Ddepp)+0.2*max(Con_p*Ddepp))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4p = Y_PG_C12_C12_4p+Y_PG_O16_C12_4p;
Y_1p = Y_PG_N14_N14_1p;
Y_6p = Y_PG_O16_O16_6p;
plot(x,Np_p/pps*Y_4p,'b'); hold on
plot(x,Np_p/pps*Y_6p,'m');
%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend('4.44 MeV', '6.13 MeV','Location', 'northwest');
[f,g]=min(Ddepp);
axis([0 (ceil(g*dx)+0.3) 0 (max(Np_p/pps*Y_4p)+0.2*max(Np_p/pps*Y_4p))]);
xlabel('Depth (cm)');
ylabel('PG/Gy/mm');
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
plot(x,Con_t*Ddept)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_t*Ddept)+0.2*max(Con_t*Ddept))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4t = Y_PG_C12_C12_4t+Y_PG_O16_C12_4t;
Y_1t = Y_PG_N14_N14_1t;
Y_2t = Y_PG_N14_N14_2t;
Y_6t = Y_PG_O16_O16_6t;
plot(x,Np_t/pps*Y_6t,'b'); hold on
plot(x,Np_t/pps*Y_4t,'m');
plot(x,Np_t/pps*Y_2t,'c');
plot(x,Np_t/pps*Y_1t,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddept);
axis([0 (ceil(g*dx)) 0 (max(Np_t/pps*Y_4t)+0.2*max(Np_t/pps*Y_4t))]);
xlabel('Depth (cm)');
ylabel('PG/Gy/mm');
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
plot(x,Con_a*Ddepa)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_a*Ddepa)+0.2*max(Con_a*Ddepa))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4a = Y_PG_C12_C12_4a+Y_PG_O16_C12_4a;
Y_1a = Y_PG_N14_N14_1a;
Y_2a = Y_PG_N14_N14_2a;
Y_6a = Y_PG_O16_O16_6a;
plot(x,Np_a/pps*Y_6a,'b'); hold on
plot(x,Np_a/pps*Y_4a,'m');
plot(x,Np_a/pps*Y_2a,'c');
plot(x,Np_a/pps*Y_1a,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepa);
axis([0 (ceil(g*dx)) 0 (max(Np_a/pps*Y_4a)+0.2*max(Np_a/pps*Y_4a))]);
xlabel('Depth (cm)');
ylabel('PG/Gy/mm');
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
plot(x,Con_b*Ddepb)
legend('Dose')
ylabel('Dose (Gy/cm³)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_b*Ddepb)+0.2*max(Con_b*Ddepb))]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_4b = Y_PG_C12_C12_4b+Y_PG_O16_C12_4b;
Y_1b = Y_PG_N14_N14_1b;
Y_2b = Y_PG_N14_N14_2b;
Y_6b = Y_PG_O16_O16_6b;
plot(x,Np_b/pps*Y_6b,'b'); hold on
plot(x,Np_b/pps*Y_4b,'m');
plot(x,Np_b/pps*Y_2b,'c');
plot(x,Np_b/pps*Y_1b,'g');

%plot(x,Y_N13p,'c')
%plot(x,Y_O15p,'m')
%plot(x,Y_C10p,'y');
legend( '6.13 MeV','4.44 MeV','2.31 MeV','1.63 MeV','Location', 'northwest');
[f,g]=min(Ddepb);
axis([0 (ceil(g*dx)) 0 (max(Np_b/pps*Y_4b)+0.2*max(Np_b/pps*Y_4b))]);
xlabel('Depth (cm)');
ylabel('PG/proton/mm');
set(gca,'FontSize',14)


%% PET activity (t) para un cierto número de protones que llegan simultaneamente
%Valores para los que se pretende calcular la actividad
calcTimes = [0 60 1000 3600]; % s

%Definimos las variables
act_C11 = zeros(numel(calcTimes), numel(x));
act_N13 = zeros(numel(calcTimes), numel(x));
act_O15 = zeros(numel(calcTimes), numel(x));
act_F18 = zeros(numel(calcTimes), numel(x));

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

deltat=1
%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % Water
    act_C11(i,:) = deltat * landa_C11 .* Y_O16_C11s .* exp(- landa_C11 * calcTimes(i));
    act_N13(i,:) = deltat * landa_N13 .* Y_O16_N13s .* exp(- landa_N13 * calcTimes(i));
    act_O15(i,:) = deltat * landa_O15 .* Y_O16_O15s .* exp(- landa_O15 * calcTimes(i));
    act_F18(i,:) = deltat * landa_F18 .* Y_O18_F18w .* exp(- landa_F18 * calcTimes(i));
    
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
    
     % Bone
    act_C11p(i,:) = deltat * landa_C11 .* Y_C11p .* exp(- landa_C11 * calcTimes(i));
    act_C10p(i,:) = deltat * landa_C10 .* Y_C10p .* exp(- landa_C10 * calcTimes(i));    
    act_N13p(i,:) = deltat * landa_N13 .* Y_N13p .* exp(- landa_N13 * calcTimes(i));
    act_O15p(i,:) = deltat * landa_O15 .* Y_O15p .* exp(- landa_O15 * calcTimes(i)); 
end
act_total = (act_C11 + act_N13 + act_O15);
act_totalt = (act_C11t + act_C10t + act_N13t + act_O15t);
act_totala = (act_C11a + act_C10a + act_N13a + act_O15a);
act_totalb = (act_C11b + act_C10b + act_N13b + act_O15b+act_Sc44b);
act_totalp = (act_C11p + act_C10p + act_N13p + act_O15p);


%% Representación gráfica de los 4 tiempos calculados
E=figure
F=figure;
G=figure;
H=figure;
I=figure;
AAA=length(calcTimes);

for i=1:numel(calcTimes)
    
        figure(E);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_w/pps*act_total(i,:),'b')
    hold on
    plot(x, Np_w/pps*act_C11(i,:),'r')
    plot(x, Np_w/pps*act_N13(i,:),'y')
    plot(x, Np_w/pps*act_O15(i,:),'c')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+/proton/mm/s ');
    legend('Total','C11','N13','O15', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
    axis([0 (ceil(g*dx)) 0 (max(Np_w/pps*act_total(i,:))+0.2*max(Np_w/pps*act_total(i,:)))]);
  
    yyaxis right
    grid on
    plot(x,Con_w*Ddep)
    ylabel('Dose (Gy/cm³)')
    set(gca,'FontSize',14)
    [f,g]=min(Ddep);
    axis([0 (ceil(g*dx)) 0 (max(Con_w*Ddep)+0.2*max(Con_w*Ddep))]);
    
    
    
    
    figure(F);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_t/pps*act_totalt(i,:),'b')
    hold on
    plot(x, Np_t/pps*act_C11t(i,:),'r')
    plot(x, Np_t/pps*act_N13t(i,:),'y')
    plot(x, Np_t/pps*act_O15t(i,:),'c')
    plot(x, Np_t/pps*act_C10t(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+/proton/mm/s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    yyaxis right
    grid on
    plot(x,Con_t*Ddept)
    ylabel('Dose (Gy/cm³)')
    set(gca,'FontSize',14)
    [f,g]=min(Ddept);
    axis([0 (ceil(g*dx)) 0 (max(Con_t*Ddept)+0.2*max(Con_t*Ddept))]);
    
    
    figure(G);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_a/pps*act_totala(i,:),'b')
    hold on
    plot(x,  Np_a/pps*act_C11a(i,:),'r')
    plot(x,  Np_a/pps*act_N13a(i,:),'y')
    plot(x,  Np_a/pps*act_O15a(i,:),'c')
    plot(x,  Np_a/pps*act_C10a(i,:),'g')
    title(sprintf('Activity at t=%i s',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+/proton/mm/s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
    axis([0 5 0 max( Np_a/pps*act_totala(i,:))+0.2*max( Np_a/pps*act_totala(i,:))]);  
    
    yyaxis right
    grid on
    plot(x,Con_a*Ddepa)
    ylabel('Dose (Gy/cm³)')
    set(gca,'FontSize',14)
    [f,g]=min(Ddept);
    axis([0 (ceil(g*dx)) 0 (max(Con_a*Ddepa)+0.2*max(Con_a*Ddepa))]);
    
    figure(H)
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_b/pps*act_totalb(i,:),'b')
    hold on
    plot(x, Np_b/pps*act_C11b(i,:),'r')
    plot(x, Np_b/pps*act_N13b(i,:),'y')
    plot(x, Np_b/pps*act_O15b(i,:),'c')
    plot(x, Np_b/pps*act_C10b(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+/Gy/mm/s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16) 
    axis([0 11 0 max(Np_b/pps*act_totalb(i,:))+0.2*max(Np_b/pps*act_totalb(i,:))]);     
    
    yyaxis right
    grid on
    plot(x,Con_b*Ddepb)
    ylabel('Dose (Gy/cm³)')
    set(gca,'FontSize',14)
    [f,g]=min(Ddepb);
    axis([0 (ceil(g*dx)) 0 (max(Con_b*Ddepb)+0.2*max(Con_b*Ddepb))]);
    
    figure(I);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_p/pps*act_totalp(i,:),'b')
    hold on
    plot(x, Np_p/pps*act_C11p(i,:),'r')
    plot(x, Np_p/pps*act_N13p(i,:),'y')
    plot(x, Np_p/pps*act_O15p(i,:),'c')
    plot(x, Np_p/pps*act_C10p(i,:),'g')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+/Gy/mm/s ');
    legend('Total','C11','N13','O15','C10', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,Con_p*Ddepp)
    ylabel('Dose (Gy/cm³)')
    set(gca,'FontSize',14)
    [f,g]=min(Ddepp);
    axis([0 (ceil(g*dx)) 0 (max(Con_p*Ddepp)+0.2*max(Con_p*Ddepp))]);
end




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
temp_totalp2=zeros(t+1);
temp_totalt2=zeros(t+1);
temp_totalb2=zeros(t+1);
temp_totala2=zeros(t+1);
temp_parcC11t=zeros(t+1,1);temp_parcC11p=zeros(t+1,1);temp_parcC11a=zeros(t+1,1);temp_parcC11b=zeros(t+1,1);
temp_parcO15t=zeros(t+1,1);temp_parcO15p=zeros(t+1,1);temp_parcO15a=zeros(t+1,1);temp_parcO15b=zeros(t+1,1);
temp_parcN13t=zeros(t+1,1);temp_parcN13p=zeros(t+1,1);temp_parcN13a=zeros(t+1,1);temp_parcN13b=zeros(t+1,1);
temp_parcC10t=zeros(t+1,1);temp_parcC10p=zeros(t+1,1);temp_parcC10a=zeros(t+1,1);temp_parcC10b=zeros(t+1,1);

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
    title('Actividad total en a lo largo del tiempo PIEL');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalt2(:,1))+0.2*max(temp_totalt2(:,1)))]); 
    
    figure;
    plot(T*deltat,(36*temp_totalp2(:,1)));
    hold on;
    plot(T*deltat,36*temp_parcO15p(:,1));
    plot(T*deltat,36*temp_parcC11p(:,1));
    plot(T*deltat,36*temp_parcN13p(:,1));
    plot(T*deltat,36*temp_parcC10p(:,1));
    %title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','O15','C11','N13','C10','Location', 'northeast');
    set(gca, 'FontSize', 16) ;
    axis([0 t*deltat 0 (max(36*temp_totalp2(:,1))+0.2*max(36*temp_totalp2(:,1)))]);
    
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
%Igual que el apartado anterior pero con otra fórmula. Es una comprobación
%de que estaba bien hecho.
cont1=0;
cont2=0;
c=1;
int_totalp2=zeros(t+1,1);
int_totalt2=zeros(t+1,1);
int_parcC11t=zeros(t+1,1);
int_parcO15t=zeros(t+1,1);
int_parcN13t=zeros(t+1,1);
int_parcC10t=zeros(t+1,1);
int_parcC11p=zeros(t+1,1);
int_parcO15p=zeros(t+1,1);
int_parcN13p=zeros(t+1,1);
int_parcC10p=zeros(t+1,1);
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
plot(x,beam_onp(1,:),'b');
plot(x,beam_offp(1,:),'r');
    title('Actividad total ');
    xlabel('z (cm)')
    ylabel('Beta+ emitters  ');
    legend('Beam On','Beam Off','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(beam_offp)+0.2*max(beam_offp))]);

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
    ylabel('Beta+ emitters ');
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
    ylabel('Beta+ emitters  ');
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
    ylabel('Beta+ emitters  ');
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
    ylabel('Beta+ emitters  ');
    legend('Beam On','Dose','Location', 'northeast');
    set(gca, 'FontSize', 16) 
    grid on;
[f,g]=min(Ddepp);
axis([ 0 (ceil(g*dx)) 0 (max(beam(1,:))+0.2*max(beam(1,:)))]);





