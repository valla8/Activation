% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es
close all

%% Par?metros a modificar:
clear all;%close all;

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control2.mat');

%PARAMETROS
dx=0.05;      %Paso del intervalo (cm)
xref=1;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=100;        %Energía inicial del haz
deltat=1;      %Inervalo de tiempo de las simulaciones
a=80/deltat;  %Tiempo de irradación del haz (s)
t=600/deltat;  %Tiempo total de la simulación
tt=240/deltat; %Tiempo de recogida de datos total
pps=3.94e10; %protones/segundo
MeVJ=1.6e-13;
%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
CsIMolecularWeight = 259.81; %g/mol
NaIMolecularWeight = 149.89; %g/mol

rho_NaI=3.67; %g/cm3
rho_CsI=4.51;

%Densidades Atomicas (A falta de las masas molares usamos
%sum(Comp_bone.*W_ele)) CORREGIR
!rho_tissue_A = 4.6243E+22; % atoms/cm3

rho_CsI_A=AvNmbr*rho_CsI/CsIMolecularWeight;
rho_NaI_A=AvNmbr*rho_NaI/NaIMolecularWeight;



%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_Cs_A = rho_CsI_A * 0.5;
rho_ICs_A = rho_CsI_A * 0.5;
rho_Na_A = rho_CsI_A * 0.5;
rho_INa_A = rho_CsI_A * 0.5;




%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
Ec = nan(size(x));
En = nan(size(x));

%Energia depositada por material
Ddepc = zeros(size(Ec));
Ddepn = zeros(size(En));

%Energia actual de cada material
currentEc = E0;
currentEn = E0;

%CREACION DE VECTORES DE YIELD
% CsI
Y_I127_Xe127c = zeros(size(x));
Y_I127_Xe125c = zeros(size(x));
Y_I127_Xe122c = zeros(size(x));
Y_I127_Xe123c = zeros(size(x));
Y_Xe127mc = zeros(size(x));

% NaI
Y_I127_Xe127n = zeros(size(x));
Y_I127_Xe125n = zeros(size(x));
Y_I127_Xe122n = zeros(size(x));
Y_I127_Xe123n = zeros(size(x));
Y_Xe127mn = zeros(size(x));
Y_Na23_Mg23n = zeros(size(x));


%CALCULO YIELD
%Recorremos cada una de las regiones del espacio caculando analíticamente
%el número de protones generado en cada una de ellas.
for i=1:(numel(x)-1)
    
    %Calculamos el poder de frenado para la energía con la que llega el
    %protón a la lámina.
    S_c = max(0,1000*S_CsI_F(currentEc*1000)); % MeV/(g/cm2)
    S_n = max(0,1000*S_NaI_F(currentEn*1000)); 
    
    %Multiplicamospor la densidad al poder de frenado
    S1c = (S_c*rho_CsI); % MeV/cm
    S1n = (S_n*rho_NaI); %*(1-Zn_fraction) + S_Znt*rho_Zn*Zn_fraction); % MeV/cm

    
    %Calculamos la energía que pierde el protón al atravesar la lámina (se
    %supone que toda la pérdida se hace al final).
    deltaEc = dx*S1c; % MeV
    deltaEn = dx*S1n; % MeV
    
    %Guardamos las energías para poder representar todo luego en función de
    %ella.
    Ec(i) = currentEc; % MeV
    En(i) = currentEn;
    
    %Calculamos la energía con la que sale de la lámina que se usará en la
    %siguiente iteración.
    currentEc = currentEc - deltaEc; % MeV
    currentEn = currentEn - deltaEn;
    
   
    
    %Guardamos la energía depositada en cada uno de los materiales.
    Ddepn(i) = deltaEn; % MeV
    Ddepc(i) = deltaEc; % MeV
    
    %Creamos unas variables con los datos anteriores para introducir a
    %continuación en los Yield.
    E1n = En(i);
    E1c = Ec(i);
    E2n = currentEn;
    E2c = currentEc;
    %Estas de abajo solo se usan por el método trapezoidal que no es el que
    %estamos teniendo en cuenta.
    E2E1n = [currentEn En(i)];
    E2E1c = [currentEc Ec(i)];
    
    %CALCULO YIELD
    %Primero se calcula la sección eficaz media en el intervalo y
    %posterioermente el yield producido. Como se introduce el número de
    %protones por segundo en realidad se calcula el yield por unidad de
    %tiempo (s).
    
    %CsI
    
    sigma_127Xe_mean = 0.5 * (max(0,I127_Xe127_F(E1c)) +max(0,I127_Xe127_F(E2c)));
    sigma_125Xe_mean = 0.5 * (max(0,I127_Xe125_F(E1c)) + max(0,I127_Xe125_F(E2c)));
    sigma_123Xe_mean = 0.5 * (max(0,I127_Xe123_F(E1c)) + max(0,I127_Xe123_F(E2c)));
    sigma_122Xe_mean = 0.5 * (max(0,I127_Xe122_F(E1c)) + max(0,I127_Xe122_F(E2c)));
    sigma_127m = 0.5 * (max(0,I127_Xem_F(E1c)) + max(0,I127_Xem_F(E2c)));
    Y_I127_Xe127c(i) =0.5 * pps * rho_CsI_A * sigma_127Xe_mean * 1e-24 * dx;
    Y_I127_Xe125c(i) =0.5 * pps * rho_CsI_A * sigma_125Xe_mean * 1e-24 * dx;
    Y_I127_Xe123c(i) =0.5 * pps * rho_CsI_A * sigma_123Xe_mean * 1e-24 * dx;
    Y_I127_Xe122c(i) =0.5 * pps * rho_CsI_A * sigma_122Xe_mean * 1e-24 * dx;
    Y_Xe127mc(i) =0.5 * pps * rho_CsI_A * sigma_127m * 1e-24 * dx; 
    
    %NaI
    
    sigma_23Na_23Mg = 0.5 * (max(0,Na23_Mg23_F(E1n)) +max(0,Na23_Mg23_F(E2n)));
    sigma_127Xe_mean = 0.5 * (max(0,I127_Xe127_F(E1n)) +max(0,I127_Xe127_F(E2n)));
    sigma_125Xe_mean = 0.5 * (max(0,I127_Xe125_F(E1n)) + max(0,I127_Xe125_F(E2n)));
    sigma_123Xe_mean = 0.5 * (max(0,I127_Xe123_F(E1n)) + max(0,I127_Xe123_F(E2n)));
    sigma_122Xe_mean = 0.5 * (max(0,I127_Xe122_F(E1n)) + max(0,I127_Xe122_F(E2n)));
    sigma_127m = 0.5 * (max(0,I127_Xem_F(E1n)) + max(0,I127_Xem_F(E2n)));
    Y_Na23_Mg23n(i) =0.5 * pps * rho_NaI_A * sigma_23Na_23Mg * 1e-24 * dx;
    Y_I127_Xe127n(i) =0.5 * pps * rho_NaI_A * sigma_127Xe_mean * 1e-24 * dx;
    Y_I127_Xe125n(i) =0.5 * pps * rho_NaI_A * sigma_125Xe_mean * 1e-24 * dx;
    Y_I127_Xe123n(i) =0.5 * pps * rho_NaI_A * sigma_123Xe_mean * 1e-24 * dx;
    Y_I127_Xe122n(i) =0.5 * pps * rho_NaI_A * sigma_122Xe_mean * 1e-24 * dx;
    Y_Xe127mn(i) =0.5 * pps * rho_NaI_A * sigma_127m * 1e-24 * dx;

    
    
end

%% Calculo de unidades Dosis

%CsI
Edepc=0;
for i=1:length(x)
    if Ddepc(i)>0.2*max(Ddepc)
        Edepc=Edepc+Ddepc(i);
    end
end
Da_c=Edepc*1.6e-13/rho_CsI;
Np_c=1/Da_c;
Con_c=Np_c*MeVJ/rho_CsI;

%Tissue
Edepn=0;
for i=1:length(x)
    if Ddepn(i)>0.2*max(Ddepn)
        Edepn=Edepn+Ddepn(i);
    end
end
Da_n=Edepn*1.6e-13/rho_NaI;
Np_n=1/Da_n;
Con_n=Np_n*MeVJ/rho_NaI;



%% Create plots de emisores beta+

x=10000*x;
dx=10000*dx;


% Figure in water
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,E)
yyaxis right
xlabel('Depth ({\mu}m)');
%hold on
plot(x,Con_c*Ddepc,'linewidth',2)
legend('Dose')
ylabel('Dose (Gy)')
set(gca,'FontSize',14)
axis([0 30 0 (max(Con_c*Ddepc)+0.2*max(Con_c*Ddepc))]);
%subplot(2,1,2)
yyaxis left
title(' CsI ');
hold on
ylabel('\beta^+ isotopes/Gy/mm³');
Y_ct = Y_I127_Xe122c+Y_I127_Xe123c+Y_I127_Xe125c+Y_I127_Xe127c+Y_Xe127mc;
plot(x,Np_c/pps*Y_I127_Xe122c,'r'); hold on
plot(x,Np_c/pps*Y_I127_Xe123c,'c')
plot(x,Np_c/pps*Y_I127_Xe125c,'m')
plot(x,Np_c/pps*Y_I127_Xe127c,'g')
plot(x,Np_c/pps*Y_Xe127mc,'b')
plot(x,Np_c/pps*Y_ct,'k','linewidth',2);
legend('122Xe','123Xe','125Xe','127Xe','127mXe','Total','Location', 'northeastoutside');
set(gca,'FontSize',14)
grid on
[f,g]=min(Ddepc);
axis([0 g*dx+0.00*g*dx  0 max(Np_c/pps*Y_ct)+0.2*max(Np_c/pps*Y_ct)]);



% Figure in NaI
figure('rend','painters','pos',[10 10 800 421])
%subplot(2,1,1)
%plot(x,E)
yyaxis right
xlabel('Depth (\mu m)');
%hold on
plot(x,Con_n*Ddepn,'linewidth',2)
legend('Dose')
ylabel('Dose (Gy)')
set(gca,'FontSize',14)
%subplot(2,1,2)
yyaxis left
title(' NaI ');
hold on
ylabel('\beta^+ isotopes/Gy/mm');
Y_nt = Y_I127_Xe122n+Y_I127_Xe123n+Y_I127_Xe125n+Y_I127_Xe127n+Y_Xe127mn+Y_Na23_Mg23n;
plot(x,Np_n/pps*Y_I127_Xe122n,'y'); hold on
plot(x,Np_n/pps*Y_I127_Xe123n,'c')
plot(x,Np_n/pps*Y_I127_Xe125n,'m')
plot(x,Np_n/pps*Y_I127_Xe127n,'g')
plot(x,Np_n/pps*Y_Xe127mn,'b')
plot(x,Np_n/pps*Y_Na23_Mg23n,'r')
plot(x,Np_n/pps*Y_nt,'k','linewidth',2);
legend('122Xe','123Xe','125Xe','127Xe','127mXe','23Mg','Total','Location', 'northeastoutside');
set(gca,'FontSize',14)
grid on
[f,g]=min(Ddepn);
axis([0 g*dx+0.00*g*dx  0 max(Np_n/pps*Y_nt)+0.2*max(Np_n/pps*Y_nt)]);





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
calcTimes = [0 10 20 60]; % s

%Definimos las variables
act_Xe127c = zeros(numel(calcTimes), numel(x));
act_Xe125c = zeros(numel(calcTimes), numel(x));
act_Xe123c = zeros(numel(calcTimes), numel(x));
act_Xe122c = zeros(numel(calcTimes), numel(x));
act_casc = zeros(numel(calcTimes), numel(x));


act_Xe127n = zeros(numel(calcTimes), numel(x));
act_Xe125n = zeros(numel(calcTimes), numel(x));
act_Xe123n = zeros(numel(calcTimes), numel(x));
act_Xe122n = zeros(numel(calcTimes), numel(x));
act_casn = zeros(numel(calcTimes), numel(x));
act_Mg23n = zeros(numel(calcTimes), numel(x));



deltat=1;
%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % CsI
    act_Xe127c(i,:) = deltat * landa_Xe127 .* Y_I127_Xe127c .* exp(- landa_Xe127 * calcTimes(i));
    act_Xe125c(i,:) = deltat * landa_Xe125 .* Y_I127_Xe125c .* exp(- landa_Xe125 * calcTimes(i));
    act_Xe123c(i,:) = deltat * landa_Xe123 .* Y_I127_Xe123c .* exp(- landa_Xe123 * calcTimes(i));
    act_Xe122c(i,:) = deltat * landa_Xe122 .* Y_I127_Xe122c .* exp(- landa_Xe122 * calcTimes(i));
    act_casc(i,:) = deltat * landa_Xe127m .* Y_Xe127mc .* exp(- landa_Xe127m * calcTimes(i));
    
    
    % CsI
    act_Xe127n(i,:) = deltat * landa_Xe127 .* Y_I127_Xe127n .* exp(- landa_Xe127 * calcTimes(i));
    act_Xe125n(i,:) = deltat * landa_Xe125 .* Y_I127_Xe125n .* exp(- landa_Xe125 * calcTimes(i));
    act_Xe123n(i,:) = deltat * landa_Xe123 .* Y_I127_Xe123n .* exp(- landa_Xe123 * calcTimes(i));
    act_Xe122n(i,:) = deltat * landa_Xe122 .* Y_I127_Xe122n .* exp(- landa_Xe122 * calcTimes(i));
    act_casn(i,:) = deltat * landa_Xe127m .* Y_Xe127mn .* exp(- landa_Xe127m * calcTimes(i));
    act_Mg23n (i,:) = deltat * landa_Mg23 .* Y_Na23_Mg23n .* exp(-landa_Mg23 *calcTimes(i));
    
end
act_totaln = (act_casn+act_Xe127n+act_Xe125n+act_Xe123n+act_Xe122n+act_Mg23n);
act_totalc = (act_casc+act_Xe127c+act_Xe125c+act_Xe123c+act_Xe122c);



%% Representación gráfica de los 4 tiempos calculados
E=figure
F=figure;
%I=figure;
AAA=length(calcTimes);

for i=1:numel(calcTimes)
    
    
    figure(E);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_c/pps*act_totalc(i,:),'k','linewidth',2)
    hold on
    plot(x, Np_c/pps*act_casc(i,:),'r-o')
    plot(x, Np_c/pps*act_Xe127c(i,:),'g-o')
    plot(x, Np_c/pps*act_Xe125c(i,:),'k-o')
    plot(x, Np_c/pps* act_Xe123c(i,:),'b-o')
    plot(x, Np_c/pps*act_Xe122c(i,:),'m-o')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth ({\mu}m)')
    ylabel('Beta+/Gy/mm³/s ');
    if i==1
    %legend('Total','Solo Agua','Xe127','N13','O15','Xe127m', 'Location', 'northwest');
    legend('Total','Xe127m', 'Xe127','Xe125','Xe123','Xe122','Location', 'northwest');
    end
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,Np_c/pps*Ddepc,'linewidth',2);
    ylabel('Dose (Gy)')
    [f,g]=min(Ddepn);
    axis([0 g*dx+0.00*g*dx 0 max(Np_c/pps*Ddepc)+0.2*max(Np_c/pps*Ddepc)]);
    
    
    
    figure(F);
    subplot(2,2,i)
    yyaxis left
    plot(x, Np_n/pps*act_totaln(i,:),'k','linewidth',2)
    hold on
    plot(x, Np_n/pps*act_casn(i,:),'r-o')
    plot(x, Np_n/pps*act_Xe127n(i,:),'g-o')
    plot(x, Np_n/pps*act_Xe125n(i,:),'k-o')
    plot(x, Np_n/pps* act_Xe123n(i,:),'y-o')
    plot(x, Np_n/pps*act_Xe122n(i,:),'m-o')
    plot(x, Np_n/pps*act_Mg23n(i,:),'b-o')
    title(sprintf('Activity at t=%i s ',calcTimes(i)));
    xlabel('Depth ({\mu}m)')
    ylabel('Beta+/Gy/mm³/s ');
    if i==1
    %legend('Total','Solo Agua','Xe127','N13','O15','Xe127m', 'Location', 'northwest');
    legend('Total','Xe127m','Xe127','Xe125','Xe123','Xe122','Mg23', 'Location', 'northwest');
    end
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,Np_n/pps*Ddepn,'linewidth',2);
    ylabel('Dose (Gy)')
    [f,g]=min(Ddepn);
    axis([0 g*dx+0.00*g*dx 0 max(Np_n/pps*Ddepn)+0.2*max(Np_n/pps*Ddepn)]);
    
    

end




%% Actividad con el tiempo
%Aquí se calcula la actividad en función del tiempo cuando el tiempo de
%irradiación es superior al delta de t, es decir, un caso clínico real
%usando un ciclotrón (no está preparado todavía para la simualción de
%pulsos de un sincrotrón). Los parámetros de tiempo de irradiación y de
%decaimiento se introducen al principio

c=1;
temp_127Xec=zeros(t+1,numel(x));temp_123Xec=zeros(t+1,numel(x));
temp_125Xec=zeros(t+1,numel(x));temp_122Xec=zeros(t+1,numel(x));
temp_127mXec=zeros(t+1,numel(x));

temp_127Xen=zeros(t+1,numel(x));temp_123Xen=zeros(t+1,numel(x));
temp_125Xen=zeros(t+1,numel(x));temp_122Xen=zeros(t+1,numel(x));
temp_127mXen=zeros(t+1,numel(x));temp_23Mgn=zeros(t+1,numel(x));

temp_parc127c=zeros(t+1,1);temp_parc125c=zeros(t+1,1);temp_parc123c=zeros(t+1,1);temp_parc122c=zeros(t+1,1);
temp_parc127mc=zeros(t+1,1);temp_totalc2=zeros(t+1,1);

temp_parc127n=zeros(t+1,1);temp_parc125n=zeros(t+1,1);temp_parc123n=zeros(t+1,1);temp_parc122n=zeros(t+1,1);
temp_parc127mn=zeros(t+1,1);temp_totaln2=zeros(t+1,1);temp_parc23=zeros(t+1,1);


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
            
            temp_127Xec(i,:)=temp_127Xec(i,:)+deltat * Y_I127_Xe127c.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1));
            temp_125Xec(i,:)=temp_125Xec(i,:)+deltat * Y_I127_Xe125c.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1));
            temp_123Xec(i,:)=temp_123Xec(i,:)+deltat * Y_I127_Xe123c.*landa_Xe123.*exp(-landa_Xe123*(deltat*j-deltat*1));
            temp_122Xec(i,:)=temp_123Xec(i,:)+deltat * Y_I127_Xe122c.*landa_Xe122.*exp(-landa_Xe122*(deltat*j-deltat*1));
            temp_127mXec(i,:)=temp_127mXec(i,:)+deltat * Y_Xe127mc.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1));
            
            temp_127Xen(i,:)=temp_127Xen(i,:)+deltat * Y_I127_Xe127n.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1));
            temp_125Xen(i,:)=temp_125Xen(i,:)+deltat * Y_I127_Xe125n.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1));
            temp_123Xen(i,:)=temp_123Xen(i,:)+deltat * Y_I127_Xe123n.*landa_Xe123.*exp(-landa_Xe123*(deltat*j-deltat*1));
            temp_122Xen(i,:)=temp_123Xen(i,:)+deltat * Y_I127_Xe122n.*landa_Xe122.*exp(-landa_Xe122*(deltat*j-deltat*1));
            temp_127mXen(i,:)=temp_127mXen(i,:)+deltat * Y_Xe127mn.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1));
            temp_23Mgn(i,:)=temp_23Mgn(i,:)+deltat * Y_Na23_Mg23n.*landa_Mg23.*exp(-landa_Mg23*(deltat*j-deltat*1));

            
           
            
        end
    else
            for j=1:a;
            temp_127Xec(i,:)=temp_127Xec(i,:)+deltat * Y_I127_Xe127c.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1+deltat*c));
            temp_125Xec(i,:)=temp_125Xec(i,:)+deltat * Y_I127_Xe125c.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1+deltat*c));
            temp_123Xec(i,:)=temp_123Xec(i,:)+deltat * Y_I127_Xe123c.*landa_Xe123.*exp(-landa_Xe123*(deltat*j-deltat*1+deltat*c));
            temp_122Xec(i,:)=temp_123Xec(i,:)+deltat * Y_I127_Xe122c.*landa_Xe122.*exp(-landa_Xe122*(deltat*j-deltat*1+deltat*c));
            temp_127mXec(i,:)=temp_127mXec(i,:)+deltat * Y_Xe127mc.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1+deltat*c));
            
            
            temp_127Xen(i,:)=temp_127Xen(i,:)+deltat * Y_I127_Xe127n.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1+deltat*c));
            temp_125Xen(i,:)=temp_125Xen(i,:)+deltat * Y_I127_Xe125n.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1+deltat*c));
            temp_123Xen(i,:)=temp_123Xen(i,:)+deltat * Y_I127_Xe123n.*landa_Xe123.*exp(-landa_Xe123*(deltat*j-deltat*1+deltat*c));
            temp_122Xen(i,:)=temp_123Xen(i,:)+deltat * Y_I127_Xe122n.*landa_Xe122.*exp(-landa_Xe122*(deltat*j-deltat*1+deltat*c));
            temp_127mXen(i,:)=temp_127mXen(i,:)+deltat * Y_Xe127mn.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1+deltat*c));
            temp_23Mgn(i,:)=temp_23Mgn(i,:)+deltat * Y_Na23_Mg23n.*landa_Mg23.*exp(-landa_Mg23*(deltat*j-deltat*1+deltat*c));
            
           
            end
            c=c+1;
    end
           
            
   % else;
           % temp_C10t(i,:)=temp_C10t(i,:)+Y_C12_C10t.*landa_C10.*exp(-landa_C10*(t-1));
            
            
            
     
    temp_totalc2(i)=sum(temp_127Xec(i,:))+sum(temp_125Xec(i,:))+sum(temp_123Xec(i,:))+sum(temp_122Xec(i,:))+sum(temp_127mXec(i,:));
    temp_parc127c(i)=sum(temp_127Xec(i,:));
    temp_parc125c(i)=sum(temp_125Xec(i,:));
    temp_parc123c(i)=sum(temp_123Xec(i,:));
    temp_parc122c(i)=sum(temp_122Xec(i,:));
    temp_parc127mc(i)=sum(temp_127mXec(i,:));            
     
    temp_totaln2(i)=sum(temp_127Xen(i,:))+sum(temp_125Xen(i,:))+sum(temp_123Xen(i,:))+sum(temp_122Xen(i,:))+sum(temp_127mXen(i,:))+sum(temp_23Mgn(i,:));
    temp_parc127n(i)=sum(temp_127Xen(i,:));
    temp_parc125n(i)=sum(temp_125Xen(i,:));
    temp_parc123n(i)=sum(temp_123Xen(i,:));
    temp_parc122n(i)=sum(temp_122Xen(i,:));
    temp_parc127mn(i)=sum(temp_127mXen(i,:));
    temp_parc23(i)=sum(temp_23Mgn(i,:));
    
   

    
    
end
    temp_totalc=temp_127Xec+temp_125Xec+temp_123Xec+temp_122Xec+temp_127mXec;
    temp_totaln=temp_127Xen+temp_125Xen+temp_123Xen+temp_122Xen+temp_127mXen+temp_23Mgn;
   
    




    %% Dibuja la Actividad total en función del tiempo.
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totalc2(:,1)),'k','linewidth',2);
    hold on;
    plot(T*deltat,temp_parc127c(:,1),'b');
    plot(T*deltat,temp_parc125c(:,1),'g');
    plot(T*deltat,temp_parc123c(:,1),'y');
    plot(T*deltat,temp_parc122c(:,1),'m');
    plot(T*deltat,temp_parc127mc(:,1),'r');
    title('Actividad total  CsI');
    grid on
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','127Xe','125Xe','123Xe','122Xe','127mXe','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totalc2(:,1))+0.2*max(temp_totalc2(:,1)))]); 
    
  
    figure;
    T=(0:t);
    plot(T*deltat,(temp_totaln2(:,1)),'k','linewidth',2);
    hold on;
    plot(T*deltat,temp_parc127n(:,1),'b');
    plot(T*deltat,temp_parc125n(:,1),'g');
    plot(T*deltat,temp_parc123n(:,1),'y');
    plot(T*deltat,temp_parc122n(:,1),'m');
    plot(T*deltat,temp_parc127mn(:,1),'r','linewidth',2);
    plot(T*deltat,temp_parc23(:,1),'c','linewidth',2);
    title('Actividad total  NaI');
    grid on
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
    legend('Actividad total','127Xe','125Xe','123Xe','122Xe','127mXe','23Mg','Location', 'northeast');
    set(gca, 'FontSize', 16); 
    axis([0 t*deltat 0 (max(temp_totaln2(:,1))+0.2*max(temp_totaln2(:,1)))]); 
      
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





