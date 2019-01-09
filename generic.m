close all
clear all

%% Par?metros a modificar:
%clear all;%close all;

% Contraste a simular
%Este programa solo representa la activación de isótopos radaictivos 
%de un contraste y los representa frente al Pico de Bragg. Para su uso se
%denomina a cada contraste con un número, y eligiendo dicho número como
%input justo aqui abajo se dibuja los resultados correspondientes. La
%asignacion de números es la siguiente.
%64Zn=1,  68Zn=2,  23Na=3, 14N=4, 127I=5, 18-O=6, 29Si=7

num=7;

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control1.mat');
load('stoppingpowers.mat');
load('CrossSections.mat');


%PARAMETROS
dx=0.01;      %Paso del intervalo (cm)
xref=15;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=50;        %Energía inicial del haz
deltat=1;      %Inervalo de tiempo de las simulaciones
a=120/deltat;  %Tiempo de irradación del haz (s)
t=900/deltat;  %Tiempo total de la simulación
tt=240/deltat; %Tiempo de recogida de datos total
pps=1.6e6; %protones/segundo
MeVJ=1.6e-13;
Con_fraction = 0.15;


%%


landa=[landa_Ga64,landa_Ga67,landa_Mg23,landa_Xe127m];
rho=[7.14 7.14 3.67 1.25 4.93 0.67*1.11 2.33 ];
AtWeight=[65.38 65.38 149.894 14 126.9 20.0276 28.05];

eners=zeros(500,500);
stop=zeros(500,500);
enerc=zeros(500,500);
cross=zeros(500,500);

eners(1:length(E_keV),1)=E_keV;
eners(1:length(E_keV),2)=E_keV;
eners(1:length(E_keV_NaI),3)=E_keV_NaI;
eners(1:length(E_keV_N),4)=E_keV_N;
eners(1:length(E_keV_Iodo),5)=E_keV_Iodo;
eners(1:length(E_keV),6)=E_keV;
eners(1:length(E_keV_Si),7)=E_keV_Si;

stop(1:length(S_Zn64),1)=S_Zn64;
stop(1:length(S_Zn68),2)=S_Zn68;
stop(1:length(S_NaI),3)=S_NaI;
stop(1:length(S_N),4)=S_N;
stop(1:length(S_Iodo),5)=S_Iodo;
stop(1:length(S_w),6)=S_w;
stop(1:length(S_Si),7)=S_Si;

cors=[length(S_Zn64) length(S_Zn68) length(S_NaI) length(S_N) length(S_Iodo) length(S_w) length(S_Si)];

enerc(1:length(Zn66_Ga66_E),1)=Zn66_Ga66_E;
enerc(1:length(Zn68_Ga68_E),2)=Zn68_Ga68_E;
enerc(1:length(Na23_Mg23_E),3)=Na23_Mg23_E;
enerc(1:length(N14_O14_E),4)=N14_O14_E;
enerc(1:length(I127_Xem_E),5)=I127_Xem_E;
enerc(1:length(O18_F18_E),6)=O18_F18_E;
enerc(1:length(Si29_P29_E),7)=Si29_P29_E;

cross(1:length(Zn66_Ga66_CS),1)=Zn66_Ga66_CS;
cross(1:length(Zn68_Ga68_CS),2)=Zn68_Ga68_CS;
cross(1:length(Na23_Mg23_CS),3)=Na23_Mg23_CS;
cross(1:length(N14_O14_CS),4)=N14_O14_CS;
cross(1:length(I127_Xem_CS),5)=I127_Xem_CS;
cross(1:length(O18_F18_CS),6)=O18_F18_CS;
cross(1:length(Si29_P29_CS),7)=Si29_P29_CS;

corc=[length(Zn66_Ga66_CS) length(Zn68_Ga68_CS) length(Na23_Mg23_CS) length(N14_O14_CS) length(I127_Xem_CS) length(O18_F18_CS) length(Si29_P29_CS)];


%Seccion Eficaz
if num==1
    Zn66_Ga66_F = fit(Zn66_Ga66_E,Zn66_Ga66_CS,'smoothingspline','SmoothingParam',0.1);
    Zn66_Ga66_F.p.coefs(1,:) = [0 0 0 0];
    Zn66_Ga66_F.p.coefs(end,:) = [0 0 0 0];
    
    Cross_F = @(x) Zn66_Ga66_F(x-3);
else
    Cross_F = fit(enerc(1:corc(num),num),cross(1:corc(num),num),'smoothingspline','SmoothingParam',0.9);
    Cross_F.p.coefs(1,:) = [0 0 0 0];
    Cross_F.p.coefs(end,:) = [0 0 0 0];
end

%Poder de Frenado
    Stop_F = fit(eners(1:cors(num),num),stop(1:cors(num),num),'smoothingspline','SmoothingParam',0.1);
    Stop_F.p.coefs(1,:) = [0 0 0 0];
    Stop_F.p.coefs(end,:) = [0 0 0 0];









%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
rho_Con=rho(num);
rho_w=1;
ConAtomicWeight=AtWeight(num);
rho_w_A =  rho_w * AvNmbr / waterMolecularWeight; % molecules / cm3
rho_Con_A = Con_fraction * rho_Con * AvNmbr / ConAtomicWeight; % molecules / cm3


%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_O16_A = rho_w_A * 0.67; % atoms/cm3


%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
E = nan(size(x));


%Energia depositada por material
Ddep = zeros(size(E));

%Energia actual de cada material
currentE = E0;

%CREACION DE VECTORES DE YIELD
% In water (full + simplified versions)
Y_Con = zeros(size(x));

%CALCULO YIELD
%Recorremos cada una de las regiones del espacio caculando analíticamente
%el número de protones generado en cada una de ellas.
for i=1:(numel(x)-1)
    
    %Calculamos el poder de frenado para la energía con la que llega el
    %protón a la lámina.
    S_w = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_Con = max(0,1000*Stop_F(currentE*1000)); % MeV/(g/cm2) 

    
    %Multiplicamospor la densidad al poder de frenado
    S_con = (S_w*rho_w)*(1-Con_fraction) + S_Con*rho_Con*Con_fraction; % MeV/cm

    
    %Calculamos la energía que pierde el protón al atravesar la lámina (se
    %supone que toda la pérdida se hace al final).
    deltaE = dx*S_con; % MeV

    
    %Guardamos las energías para poder representar todo luego en función de
    %ella.
    E(i) = currentE; % MeV

    
    %Calculamos la energía con la que sale de la lámina que se usará en la
    %siguiente iteración.
    currentE = currentE - deltaE; % MeV

   
    %Guardamos la energía depositada en cada uno de los materiales.
    Ddep(i) = deltaE; % MeV

    
    %Creamos unas variables con los datos anteriores para introducir a
    %continuación en los Yield.
    E1 = E(i);
    E2 = currentE;

    %Estas de abajo solo se usan por el método trapezoidal que no es el que
    %estamos teniendo en cuenta.
    E2E1 = [currentE E(i)];

    
    %CALCULO YIELD
    %Primero se calcula la sección eficaz media en el intervalo y
    %posterioermente el yield producido. Como se introduce el número de
    %protones por segundo en realidad se calcula el yield por unidad de
    %tiempo (s).
    
    % Water (full + simplified)
    sigma_Con =  (max(0,Cross_F(E1)) + max(0,Cross_F(E2)));
    Y_Con(i) = rho_Con_A * sigma_Con * 1e-24 * dx;
    
  
    
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
ylabel('Dose (Gy*cm³)')
set(gca,'FontSize',14)
%subplot(2,1,2)
yyaxis left
title(' Water ');
hold on
ylabel('\beta^+ isotopes/Gy/mm');
plot(x,Np_w/pps*Y_Con,'k','linewidth',2);
legend('Contraste','Total','Location', 'northeastoutside');
set(gca,'FontSize',14)
grid on
[f,g]=min(Ddep);
axis([g*dx-1 g*dx  0 max(Np_w/pps*Y_Con)+0.2*max(Np_w/pps*Y_Con)]);







