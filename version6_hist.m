% C?digo de prueba de c?lculo de actividad con Zn
% (C)Victor Valladolid Onecha 2019
% vicvalla@ucm.es
close all

%% Par?metros a modificar:
clear all;%close all;

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control2.mat');
landa_F18 =  log(2) / 6586;

%PARAMETROS
dx=0.05;      %Paso del intervalo (cm)
xref=10;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=100;        %Energía inicial del haz
deltat=1;      %Inervalo de tiempo de las simulaciones
a=80/deltat;  %Tiempo de irradación del haz (s)
t=600/deltat;  %Tiempo total de la simulación
tt=240/deltat; %Tiempo de recogida de datos total
pps=1e6; %protones/segundo
MeVJ=1.6e-13;
O18_fraction=0.15;

%% Densidades Atómicas

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
PMMA_Molar=100.12; %g/mol
rho_w = 1; % g/cm3
rho_w18 = 1.1; % g/cm3

W_ele=[1.0079 12.0110 14.0067 15.994 18];

Comp_water = [0.667 0 0 0.333 0];
Comp_h2o18 = [0.667 0 0 0 0.333];
Comp_water_h2o18 = O18_fraction*Comp_h2o18 + (1-O18_fraction)* Comp_water;


%Densidades Atomicas 

!rho_tissue_A = 4.6243E+22; % atoms/cm3
%rho_PMMA_A = AvNmbr*rho_PMMA/PMMA_Molar;  % atoms/cm3
%rho_w_A =  rho_8w * AvNmbr / waterMolecularWeight; % molecules / cm3
rho_w_A =  (1-O18_fraction) * rho_w * AvNmbr / sum(Comp_water.*W_ele); % molecules / cm3
rho_w18_A =  O18_fraction * rho_w18 * AvNmbr / sum(Comp_water_h2o18.*W_ele); % molecules / cm3

%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_O16_A = rho_w_A * Comp_water(4) * O16_ab; % atoms/cm3
rho_O18_A = rho_w_A * Comp_water_h2o18(5); % atoms/cm3



%% Calculo Analítico

%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
E = nan(size(x));
Et = nan(size(x));
Eb = nan(size(x));
Ea = nan(size(x));
Ep = nan(size(x));




%CREACION DE VECTORES DE YIELD
% In water (full + simplified versions)
Y_O16_C11s = zeros(size(x));
Y_O16_N13s = zeros(size(x));
Y_O16_O15s = zeros(size(x));
Y_O18_F18w = zeros(size(x));
% Y_PG_C12_C12_4w = zeros(size(x));
% Y_PG_O16_C12_4w = zeros(size(x));
% Y_PG_N14_N14_1w = zeros(size(x));
% Y_PG_N14_N14_2w = zeros(size(x));
% Y_PG_O16_O16_6w = zeros(size(x));




%CALCULO YIELD
%Recorremos cada una de las regiones del espacio caculando analíticamente
%el número de protones generado en cada una de ellas.

%Metemos el histograma

     row=2000;  col=201;
     fin=fopen('Tot.raw','r');
     I=fread(fin,row*col,'single'); 
     histo=reshape(I,row,col);
     histo=histo/4.6606e+05;

for i=2:201
        ii=i-1
    

    %CALCULO YIELD
    % En esta versión introducimos los histogramas de energía de cada bin
    % obtenidos es TOPAS, por lo que para cada bin tenemos que meter un
    % bucle de los 2000 bines de cada histograma
    for j=10:2000
    jj=j/10;
    
    % Water (full + simplified)
%     sigma_PG_C12_4w = max(0,PG_C12_C12_4_F(E(j)))
%     sigma_PG_O16_4w = max(0,PG_O16_C12_4_F(E1)) + max(0,PG_O16_C12_4_F(E2)));
%     sigma_PG_N14_1w = max(0,PG_N14_N14_1_F(E1)) + max(0,PG_N14_N14_1_F(E2)));
%     sigma_PG_N14_2w = max(0,PG_N14_N14_2_F(E1)) + max(0,PG_N14_N14_2_F(E2)));
%     sigma_PG_O16_6w = max(0,PG_O16_O16_6_F(E1)) + max(0,PG_O16_O16_6_F(E2)));
    sigma_C11_mean = histo(j,i) * max(0,O16_C11_F(jj));
    sigma_N13_mean = histo(j,i) * max(0,O16_N13_F(jj));
    sigma_O15_mean = histo(j,i) * max(0,O16_O15_F(jj));
    sigma_F18_mean = histo(j,i) * max(0,O18_F18_F(jj));
    Y_O16_C11s(ii) = Y_O16_C11s(ii) + rho_O16_A * sigma_C11_mean * 1e-24 * dx;
    Y_O16_N13s(ii) = Y_O16_N13s(ii) + rho_O16_A * sigma_N13_mean * 1e-24 * dx;
    Y_O16_O15s(ii) = Y_O16_O15s(ii) + rho_O16_A * sigma_O15_mean * 1e-24 * dx;
    Y_O18_F18w(ii) = Y_O18_F18w(ii) + rho_O18_A * sigma_F18_mean * 1e-24 * dx;
%     Y_PG_O16_C12_4w(i) = pps * rho_O16_A * sigma_PG_O16_4w * 1e-24 * dx;
%     Y_PG_O16_O16_6w(i) = pps * rho_O16_A * sigma_PG_O16_6w * 1e-24 * dx;
  
    end
    
    
end



%%

AA=zeros(200,5);
AA(:,1)=linspace(1,200,200);
AA(:,2)=Y_O16_O15s(1:200);
AA(:,3)=Y_O16_N13s(1:200)/1000;
AA(:,4)=Y_O16_C11s(1:200);
AA(:,5)=1000*Y_O18_F18w(2:201);
AA(:,6)=Y_O16_O15s(1:200)+Y_O16_N13s(1:200)/1000+Y_O16_C11s(1:200)+1000 *Y_O18_F18w(1:200);

%%
plot(AA(:,1),Y_O16_O15s(2:201));
hold on;
plot(AA(:,1),Y_O16_N13s(2:201)/1000);
plot(AA(:,1),Y_O16_C11s(2:201));
plot(AA(:,1),1000*Y_O18_F18w(2:201));
%plot(AA(:,1),

