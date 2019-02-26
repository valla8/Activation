% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es
close all
clear all;%close all;

%% Par?metros a modificar:

%Cargamos las vidas medias, secciones eficaces y stopping power para
%ahorrar tiempo de calculo.
load('control2.mat');

%PARAMETROS
dx=0.1;      %Paso del intervalo (cm)
xref=10;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=100;        %Energía inicial del haz

pps=6.25e9; %protones/segundo
MeVJ=1.6e-13;
landa_F18 =  log(2) / 6586;
O18_fraction=0.0;
%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
PMMA_Molar=100.12; %g/mol
rho_w = 1; % g/cm3
rho_w18 = 1.11; % g/cm3
rho_Zn = 7.14; %g/cm3
rho_tissue = 1.1; %g/cm3
rho_bo = 1.85; %g/cm3
rho_adipose = 0.92; %g/cm3
rho_PMMA= 1.18; %g/cm3
rho_bet =1.025 ; %g/cm3
rho_CsI=4.51;

W_ele=[1.0079 12.0110 14.0067 15.994 18];

Comp_water = [0.667 0 0 0.333  0];
Comp_h2o18 = [0.667 0 0 0 0.333];
Comp_water_h2o18 = O18_fraction*Comp_h2o18 + (1-O18_fraction)* Comp_water;

rho_w_A =  (1-O18_fraction) * rho_w * AvNmbr / sum(Comp_water.*W_ele); % molecules / cm3
rho_w18_A =  O18_fraction * rho_w18 * AvNmbr / sum(Comp_water_h2o18.*W_ele); % molecules / cm3


%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_O16_A = rho_w_A * Comp_water_h2o18(4) * O16_ab; % atoms/cm3
rho_O18_A = rho_w18_A * Comp_water_h2o18(5); % atoms/cm3

%%

%Se generan los vectores que se usarán en la simulación
x = 0:dx:xref; % posiciones en cm.
E = zeros(size(x));

%Energia depositada por material
Ddep = zeros(size(E));

%Energia actual de cada material
currentE = E0;

%CREACION DE VECTORES DE YIELD
% In water (full + simplified versions)


%CALCULO YIELD
%Recorremos cada una de las regiones del espacio caculando analíticamente
%el número de protones generado en cada una de ellas.
for i=1:(numel(x)-1)
    
    %Calculamos el poder de frenado para la energía con la que llega el
    %protón a la lámina.
    S_w = max(0,1000*S_w_F(currentE*1000)); % MeV/(g/cm2)
    S_w18 = max(0,1000*S_w18_F(currentE*1000)); % MeV/(g/cm2)
    
    %Multiplicamospor la densidad al poder de frenado
    S1 = (1-O18_fraction)*S_w*rho_w+O18_fraction*S_w18*rho_w18; % MeV/cm
    
    %Calculamos la energía que pierde el protón al atravesar la lámina (se
    %supone que toda la pérdida se hace al final).
    deltaE = dx*S1; % MeV
    
    %Guardamos las energías para poder representar todo luego en función de
    %ella.
    E(i) = currentE; % MeV
    
    if E(i) < 0
        E(i)=0;
    end
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
    
end
        
      %%
    
    %CALCULO YIELD
    %Primero se calcula la sección eficaz media en el intervalo y
    %posterioermente el yield producido. Como se introduce el número de
    %protones por segundo en realidad se calcula el yield por unidad de
    %tiempo (s).
    Y_O16_C11s = zeros(size(x));
    Y_O16_N13s = zeros(size(x));
    Y_O16_O15s = zeros(size(x));
    Y_O18_F18w = zeros(size(x));

Fl=zeros(size(x));
Fl2=zeros(size(x));

    EE=zeros(size(x));
        
    [f,g]=min(E);
    R0=g*dx
    
    p=1.5;
    alpha=1.6*10^(-3);
    sigee=0.7
    
    sige=sigee*alpha^(1/p)*p*R0^(1-1/p);
    sigr=0.012*(alpha*E0^p)^0.935;
    sig2=sige^2+sigr^2;
    sig=sqrt(sig2);
    
     for i=1:length(x)
   
        for k=1:length(x)

         z=dx*i;
         beta=0.012; %cm-1
         resr=R0-dx*k;
         if resr >= 0 
 %        Fl(i)=exp(-(resr-(R0-z))^2/(2*sig* 2));    
 %        Fl(k)=1/(19.9706*sqrt(2*3.1416)*sig)*(1+beta*resr)/(1+beta*R0)*exp(-(resr-(R0-z))^2/(2*sig2))
         Fl(k)=1/(11.9809*sqrt(2*3.1416)*sig)*(1+beta*resr)/(1+beta*R0)*exp(-(resr-(R0-z))^2/(2*sig2));
 %        1/(11.9809*sqrt(2*3.1416)*sig)*(1+beta*resr)/(1+beta*R0)*exp(-(resr-(R0-z))^2/(2*sig2))
 %       Fl(k)=(1+beta*resr)/(1+beta*R0)*exp(-(resr-(R0-z))^2/(2*sig2));
         else
         Fl(k)=0; 
         end
         
       
         Fl2(i)=Fl2(i)+Fl(k);
         
         
    
    % Water (full + simplified)
    sigma_C11_mean = max(0,O16_C11_F(E(k)));
    sigma_N13_mean = max(0,O16_N13_F(E(k))) ;
    sigma_O15_mean = max(0,O16_O15_F(E(k)));
    sigma_F18_mean = max(0,O18_F18_F(E(k)));
    Y_O16_C11s(i) = Y_O16_C11s(i) +  Fl(k)*rho_O16_A * sigma_C11_mean * 1e-24 * dx;
    Y_O16_N13s(i) = Y_O16_N13s(i) +  Fl(k)*rho_O16_A * sigma_N13_mean * 1e-24 * dx;
    Y_O16_O15s(i) = Y_O16_O15s(i) +  Fl(k)*rho_O16_A * sigma_O15_mean * 1e-24 * dx;
%     k;
%     Fl(k)
%     rho_O16_A * sigma_O15_mean * 1e-24 * dx
    
    Y_O18_F18w(i) = Y_O18_F18w(i) +  Fl(k)*rho_O18_A * sigma_F18_mean * 1e-24 * dx;
    
%    EE(i) = EE(i) ;
    
    
    
        end

    
    
end
    Y_N13=Y_O16_N13s;
    Y_F18=Y_O18_F18w;
    Y_O15=Y_O16_O15s;
%%




figure
hold on
plot(E,O16_O15_F(E))
plot(E,O16_N13_F(E))
plot(E,O16_C11_F(E))
plot(E,O18_F18_F(E))

legend('15-O','13-N','11-C','18-F')




%%
figure
hold on

% plot(x,Ddep/1e6,'k')
 plot(x,Y_N13)
 plot(x,Y_O15)
 plot(x,Y_F18)

legend('Dose','13-N','18-F')

 %% PET activity (t) para un cierto número de protones que llegan simultaneamente
% %Valores para los que se pretende calcular la actividad
calcTimes = [0 30 60 120 3600*2 ]; % s

%Definimos las variables
act_13N = zeros(numel(calcTimes), numel(x));
act_18F = zeros(numel(calcTimes), numel(x));


nC=1/1.6e-10;



deltat=1
%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % Water
    act_13N(i,:) = deltat * landa_N13 .* Y_N13 .* exp(- landa_N13 * calcTimes(i));
    act_18F(i,:) = deltat * landa_F18 .* Y_F18 .* exp(- landa_F18 * calcTimes(i));
    

   
end
act_total = (act_13N+act_18F);
%act_totalp = (act_C11p + act_C10p + act_N13p + act_O15p);
yyaxis left
plot(x,act_total(1,:),'k')
hold on
plot(x,act_total(2,:),'r')
plot(x,act_total(3,:),'b')
plot(x,act_total(4,:),'y')
plot(x,act_total(5,:),'g')
grid on
[f,g]=min(Ddep);
axis([0 g*dx  0 (max(act_total(1,:))+0.2*max(act_total(1,:)))]);
yyaxis right
plot(x,Ddep)





%%
%Definimos las variables
act_13Nt = zeros(1,360);
act_18Ft = zeros(1,360);
TT=[1:360];

for i=1:360


    act_13Nt(i) =  landa_N13 .* sum(Y_F18) .* exp(- landa_N13 * i);
    act_18Ft(i) =  landa_F18 .* sum(Y_N13) .* exp(- landa_F18 * i);
    
end


plot(TT,act_13Nt)
hold on
plot(TT,act_18Ft)
legend('13-N','18-F')


%% Actividad con el tiempo
%Aquí se calcula la actividad en función del tiempo cuando el tiempo de
%irradiación es superior al delta de t, es decir, un caso clínico real
%usando un ciclotrón (no está preparado todavía para la simualción de
%pulsos de un sincrotrón). Los parámetros de tiempo de irradiación y de
%decaimiento se introducen al principio
deltat=1;      %Inervalo de tiempo de las simulaciones
a=1800/deltat;  %Tiempo de irradación del haz (s)
t=4000/deltat;  %Tiempo total de la simulación


c=1;
temp_13N=zeros(t+1,numel(x));
temp_18F=zeros(t+1,numel(x));

temp_total2=zeros(t+1,1);
temp_parc13N=zeros(t+1,1);
temp_parc18F=zeros(t+1,1);

%CALCULO
%Para hacer el calculo suponemos que cada intervalod de tiempo llegan un
%número de protones que se frenan inmediatamente y genera los isótopos
%correspondientes. En cada iteración se calcula el número de los isótopos
%de las iteraciones anteriores que se ha desintegrado y además si estamos
%en el tiempo de irradiación se suma la contribución de los protones que
%llegan.
d=0
buff=0
f=0;

for i=1:t+1
    b=i;
%    if buff<buff1
%    if i<tau*2.5/deltat
    if i<a
           d=d+1;
        for j=1:b
            
            temp_18F(i,:)=temp_18F(i,:)+pps*deltat * Y_F18.*landa_F18.*exp(-landa_F18*(deltat*j-deltat*1));
            temp_13N(i,:)=temp_13N(i,:)+pps*deltat * Y_N13.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1));

            
        end
    else
            for j=1:d;
                
             temp_18F(i,:)=temp_18F(i,:)+pps*deltat *Y_F18.*landa_F18.*exp(-landa_F18*(deltat*j-deltat*1+deltat*c));
             temp_13N(i,:)=temp_13N(i,:)+pps*deltat *Y_N13.*landa_N13.*exp(-landa_N13*(deltat*j-deltat*1+deltat*c));


            end
            c=c+1;

            
  %          sum(temp_Ba133m(i,:))
    end
           
            
   % else;
           % temp_C10t(i,:)=temp_C10t(i,:)+Y_C12_C10t.*landa_C10.*exp(-landa_C10*(t-1));
            
            
            
     
    temp_total2(i)=sum(sum(temp_13N(i,:))+sum(temp_18F(i,:)));
    temp_parc18F(i)=sum(temp_18F(i,:));
    temp_parc13N(i)=sum(temp_13N(i,:));


    
    
end
    temp_total=temp_18F+temp_13N;
    




    %% Dibuja la Actividad total en función del tiempo.
    figure;
    set(gca, 'FontSize', 16); 
    T=(0:t);
    hold on;
%    tottot=0.003122*0.38*temp_parcXem127+0.00285*0.68*temp_parcXem127+0.002535*0.1769*temp_parcBa133m;
%    plot(T*deltat,tottot/1000,'k','linewidth',2);
%    plot(T*deltat,(0.003122*0.38*temp_parcXem127)/1000,'r','linewidth',2);
%    plot(T*deltat,(0.00285*0.68*temp_parcXem127)/1000,'g','linewidth',2);
    plot(T*deltat,(temp_parc18F)/1000,'r');
    plot(T*deltat,(temp_parc13N)/1000,'g');
    grid on
    title('Actividad/nC ');
    xlabel('Tiempo (s)');
    ylabel('Actividad (kBq)');
    legend('18-F','13-N');
 %   legend('Xe127*-172 keV','Xe127*-124 keV','Xe127','Xe125','Xe125','Xe123','Xe122');
%    axis([0 t*deltat 0 (0.68*max(temp_total2(:,1)/1000)+0.2*max(temp_total2(:,1)/1000))]); 
%     
max(temp_parc18F)/1000
max(temp_parc13N)/1000
