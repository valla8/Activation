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
dx=0.001;      %Paso del intervalo (cm)
xref=4.0;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=100;        %Energía inicial del haz

pps=6.25e9; %protones/segundo
MeVJ=1.6e-13;
landa_F18 =  log(2) / 6586;
landa_Ba133 = log(2) / 3.327e8;
landa_Ba133m = log(2) / 140148;
O18_fraction=0.00;
%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
waterMolecularWeight = 18.01528; %g/mol
PMMA_Molar=100.12; %g/mol
rho_w = 1; % g/cm3
rho_w18 = 1.1; % g/cm3
rho_Zn = 7.14; %g/cm3
rho_tissue = 1.1; %g/cm3
rho_bone = 1.85; %g/cm3
rho_adipose = 0.92; %g/cm3
rho_PMMA= 1.18; %g/cm3
rho_bet =1.025 ; %g/cm3
rho_CsI=4.51;

W_ele=[ 126.9 132.9 ];

Comp_CsI=[0.5 0.5]

%Densidades Atomicas (A falta de las masas molares usamos
%sum(Comp_bone.*W_ele)) CORREGIR
% !rho_tissue_A = 4.6243E+22; % atoms/cm3
% rho_tissue_A = AvNmbr*rho_tissue/sum(Comp_tissue.*W_ele);  % atoms/cm3
% rho_bone_A = AvNmbr*rho_bone/sum(Comp_bone.*W_ele);  % atoms/cm3
% rho_adipose_A = AvNmbr*rho_adipose/sum(Comp_adipose.*W_ele);  % atoms/cm3
% rho_PMMA_A = AvNmbr*rho_PMMA/sum(Comp_PMMA.*W_ele);  % atoms/cm3
%rho_PMMA_A = AvNmbr*rho_PMMA/PMMA_Molar;  % atoms/cm3
%rho_w_A =  rho_w * AvNmbr / waterMolecularWeight; % molecules / cm3
%rho_w_A =  (1-O18_fraction) * rho_w * AvNmbr / ; % molecules / cm3
%rho_w18_A =  O18_fraction * rho_w18 * AvNmbr / sum(Comp_water_h2o18.*W_ele); % molecules / cm3
% rho_w18_A =  O18_fraction * rho_w18 * AvNmbr / sum(Comp_water_h2o18.*W_ele); % molecules / cm3
% ZnAtomicWeight = 65.38; % g/mol
% rho_Zn_A = Zn_fraction * rho_Zn * AvNmbr / ZnAtomicWeight; % molecules / cm3
rho_CsI_A =  rho_CsI * AvNmbr / sum(Comp_CsI.*W_ele); % molecules / cm3

%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_Cs_A = rho_CsI_A * Comp_CsI(1); % atoms/cm3
rho_I_A = rho_CsI_A * Comp_CsI(2); % atoms/cm3

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
    S_CsI = max(0,1000*S_CsI_F(currentE*1000)); % MeV/(g/cm2)
    
    %Multiplicamospor la densidad al poder de frenado
    S1 = S_CsI*rho_CsI; % MeV/cm
    
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
    Y_127I_127Xe = zeros(size(x));
    Y_127I_125Xe = zeros(size(x));
    Y_127I_122Xe = zeros(size(x));
    Y_127I_123Xe = zeros(size(x));
    Y_127Xem = zeros(size(x));
    Y_133Cs_133Ba = zeros(size(x));
    Y_133Cs_133Bam = zeros(size(x));
    
    CCC=0;
    
     for i=1:length(x)-1

     
         
    
    % Water (full + simplified)

   

%     Y_127Xem(i) = Y_127Xem(i)  +  rho_I_A * max(0,(I127_Xem_F(E(i))+I127_Xem_F(E(i+1)))/2) * 1e-24 * dx;
     Y_127Xem(i) = Y_127Xem(i)  +  rho_I_A * max(0,1e3*(I127_Xe127_F(E(i))+I127_Xe127_F(E(i+1)))/4) * 1e-24 * dx;
     Y_127I_127Xe(i) = Y_127I_127Xe(i)  +  rho_I_A * max(0,1000*I127_Xe127_F(E(i))) * 1e-24 * dx;
        if i<length(x)
        CCC=CCC+(-E(i+1)+E(i))*max(0,I127_Xem_F(E(i)));
        end
     Y_127I_125Xe(i) = Y_127I_125Xe(i)  +  rho_I_A * max(0,1e3*I127_Xe125_F(E(i))) * 1e-24 * dx;
     Y_127I_123Xe(i) = Y_127I_123Xe(i)  +  rho_I_A * max(0,1e3*I127_Xe123_F(E(i))) * 1e-24 * dx;
     Y_127I_122Xe(i) = Y_127I_122Xe(i)  +  rho_I_A * max(0,1e3*I127_Xe122_F(E(i))) * 1e-24 * dx;
     Y_133Cs_133Ba(i) = Y_133Cs_133Ba(i)  +  rho_Cs_A * max(0,Cs133_Ba133_F(E(i))) * 1e-24 * dx;
     Y_133Cs_133Bam(i) = Y_133Cs_133Bam(i)  +  rho_Cs_A * max(0,Cs133_Ba133m_F(E(i))) * 1e-24 * dx;
    
%     EE(i) = EE(i) + Fl(k) * Ddep(k);
    
    
    
     end
        
     

    Y_Xe127=Y_127I_127Xe;
    Y_Xe125=Y_127I_125Xe;
    Y_Xe123=Y_127I_123Xe;
    Y_Xe122=Y_127I_122Xe;
    Y_Ba133 = Y_133Cs_133Ba;
    Y_Ba133m = Y_133Cs_133Bam;

    
    
    % end
%%




figure
hold on
plot(E,I127_Xem_F(E))
plot(E,1e3*I127_Xe122_F(E))
plot(E,1e3*I127_Xe123_F(E))
plot(E,1e3*I127_Xe125_F(E))
plot(E,10e2*I127_Xe127_F(E))
plot(E,Cs133_Ba133_F(E))
plot(E,Cs133_Ba133m_F(E))
legend('127Xem','122Xe','123Xe','125Xe','127Xe','133Ba','133mBa')




%%
figure
hold on

plot(x,Ddep/1e5,'k')
 plot(x,Y_127Xem)
 plot(x,Y_127I_122Xe)
 plot(x,Y_127I_123Xe)
 plot(x,Y_127I_125Xe)
 plot(x,Y_127I_127Xe)
 plot(x,Y_133Cs_133Ba)
 plot(x,Y_133Cs_133Bam)
legend('Dose','127Xem','122Xe','123Xe','125Xe','127Xe','133Ba','133mBa')

 %% PET activity (t) para un cierto número de protones que llegan simultaneamente
% %Valores para los que se pretende calcular la actividad
calcTimes = [0 30 60 120 3600 ]; % s

%Definimos las variables
act_Xe127m = zeros(numel(calcTimes), numel(x));
act_Xe127 = zeros(numel(calcTimes), numel(x));
act_Xe125 = zeros(numel(calcTimes), numel(x));
act_Xe123 = zeros(numel(calcTimes), numel(x));
act_Xe122 = zeros(numel(calcTimes), numel(x));
act_Ba133 = zeros(numel(calcTimes), numel(x));
act_Ba133m = zeros(numel(calcTimes), numel(x));


nC=1/1.6e-10;



deltat=1
%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % Water
    act_Xe127m(i,:) = deltat * landa_Xe127m .* Y_127Xem .* exp(- landa_Xe127m * calcTimes(i));
    act_Xe127(i,:) = deltat * landa_Xe127 .* Y_Xe127 .* exp(- landa_Xe127 * calcTimes(i));
    act_Xe125(i,:) = deltat * landa_Xe125 .* Y_Xe125 .* exp(- landa_Xe125 * calcTimes(i));
    act_Xe123(i,:) = deltat * landa_Xe123 .* Y_Xe123 .* exp(- landa_Xe123 * calcTimes(i));
    act_Xe122(i,:) = deltat * landa_Xe122 .* Y_Xe122 .* exp(- landa_Xe122 * calcTimes(i));
    act_Ba133(i,:) = deltat * landa_Ba133 .* Y_Ba133 .* exp(- landa_Ba133 * calcTimes(i));
    act_Ba133m(i,:) = deltat * landa_Ba133m .* Y_Ba133m .* exp(- landa_Ba133m * calcTimes(i));
    

   
end
act_total = (act_Xe127m+act_Xe125+act_Xe123+act_Xe122+act_Ba133);
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


BB=zeros(length(x),length(calcTimes)+1);
BB(:,1)=linspace(0,max(x),length(x));
BB(:,2)=nC*act_total(1,:);
BB(:,3)=nC*act_total(2,:);
BB(:,4)=nC*act_total(3,:);
BB(:,5)=nC*act_total(4,:);
BB(:,6)=nC*act_total(5,:);


%%
%Definimos las variables
act_Xe127mt = zeros(1,360);
act_Xe127t = zeros(1,360);
act_Xe125t = zeros(1,360);
act_Xe123t = zeros(1,360);
act_Xe122t = zeros(1,360);
act_Ba133t = zeros(1,360);
act_Ba133mt = zeros(1,360);
TT=[1:360];

for i=1:360


    act_Xe127mt(i) =  landa_Xe127m .* sum(Y_127Xem) .* exp(- landa_Xe127m * i);
    act_Xe127t(i) =  landa_Xe127 .* sum(Y_Xe127) .* exp(- landa_Xe127m * i);
    act_Xe125t(i) =  landa_Xe125 .* sum(Y_Xe125) .* exp(- landa_Xe125 * i);
    act_Xe123t(i) =  landa_Xe123 .* sum(Y_Xe123) .* exp(- landa_Xe123 * i);
    act_Xe122t(i) =  landa_Xe122 .* sum(Y_Xe122) .* exp(- landa_Xe123 * i);
    act_Ba133t(i) =  landa_Ba133 .* sum(Y_Ba133) .* exp(- landa_Ba133 * i);
    act_Ba133mt(i) =  landa_Ba133m .* sum(Y_Ba133m) .* exp(- landa_Ba133m * i);
    
end


plot(TT,act_Xe127mt)
hold on
plot(TT,act_Xe127t)
plot(TT,act_Xe125t)
plot(TT,act_Xe123t)
plot(TT,act_Xe122t)
plot(TT,act_Ba133t)
plot(TT,act_Ba133mt)
legend('127Xem','127Xe','Ba133','Ba133m')


%% Actividad con el tiempo
%Aquí se calcula la actividad en función del tiempo cuando el tiempo de
%irradiación es superior al delta de t, es decir, un caso clínico real
%usando un ciclotrón (no está preparado todavía para la simualción de
%pulsos de un sincrotrón). Los parámetros de tiempo de irradiación y de
%decaimiento se introducen al principio
deltat=1;      %Inervalo de tiempo de las simulaciones
a=600/deltat;  %Tiempo de irradación del haz (s)
t=600/deltat;  %Tiempo total de la simulación


c=1;
temp_Xem127=zeros(t+1,numel(x));
temp_Xe127=zeros(t+1,numel(x));
temp_Xe125=zeros(t+1,numel(x));
temp_Xe123=zeros(t+1,numel(x));
temp_Xe122=zeros(t+1,numel(x));
temp_Ba133=zeros(t+1,numel(x));
temp_Ba133m=zeros(t+1,numel(x));

temp_total2=zeros(t+1,1);
temp_parcXem127=zeros(t+1,1);
temp_parcXe127=zeros(t+1,1);
temp_parcXe125=zeros(t+1,1);
temp_parcXe123=zeros(t+1,1);
temp_parcXe122=zeros(t+1,1);
temp_parcBa133=zeros(t+1,1);
temp_parcBa133m=zeros(t+1,1);

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
buff1=0.9*sum(Y_127Xem);
tau=log(2)/landa_Xe127m
for i=1:t+1
    b=i;
%    if buff<buff1
    if i<tau*2.5/deltat
           d=d+1;
    %if i<a
        for j=1:b
            
            temp_Xem127(i,:)=temp_Xem127(i,:)+pps*deltat * Y_127Xem.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1));
            temp_Xe127(i,:)=temp_Xe127(i,:)+pps*deltat * Y_Xe127.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1));
            temp_Xe125(i,:)=temp_Xe125(i,:)+pps*deltat * Y_Xe125.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1));
            temp_Xe123(i,:)=temp_Xe123(i,:)+pps*deltat * Y_Xe123.*landa_Xe123.*exp(-landa_Xe125*(deltat*j-deltat*1));
            temp_Xe122(i,:)=temp_Xe122(i,:)+pps*deltat * Y_Xe122.*landa_Xe122.*exp(-landa_Xe125*(deltat*j-deltat*1));
            temp_Ba133(i,:)=temp_Ba133(i,:)+pps*deltat * Y_Ba133.*landa_Ba133.*exp(-landa_Ba133*(deltat*j-deltat*1));
            temp_Ba133m(i,:)=temp_Ba133m(i,:)+pps*deltat * Y_Ba133m.*landa_Ba133m.*exp(-landa_Ba133m*(deltat*j-deltat*1));
           buff=sum(temp_Xem127(i,:))/pps;
            
        end
    else
            for j=1:d;
                
            temp_Xem127(i,:)=temp_Xem127(i,:)+pps*deltat *Y_127Xem.*landa_Xe127m.*exp(-landa_Xe127m*(deltat*j-deltat*1+deltat*c));
             temp_Xe127(i,:)=temp_Xe127(i,:)+pps*deltat *Y_Xe127.*landa_Xe127.*exp(-landa_Xe127*(deltat*j-deltat*1+deltat*c));
            temp_Xe125(i,:)=temp_Xe125(i,:)+pps*deltat *Y_Xe125.*landa_Xe125.*exp(-landa_Xe125*(deltat*j-deltat*1+deltat*c));
            temp_Xe123(i,:)=temp_Xe123(i,:)+pps*deltat *Y_Xe123.*landa_Xe123.*exp(-landa_Xe123*(deltat*j-deltat*1+deltat*c));
            temp_Xe122(i,:)=temp_Xe122(i,:)+pps*deltat *Y_Xe122.*landa_Xe122.*exp(-landa_Xe122*(deltat*j-deltat*1+deltat*c));
            temp_Ba133(i,:)=temp_Ba133(i,:)+pps*deltat *Y_Ba133.*landa_Ba133.*exp(-landa_Ba133*(deltat*j-deltat*1+deltat*c));
            temp_Ba133m(i,:)=temp_Ba133m(i,:)+pps*deltat *Y_Ba133m.*landa_Ba133m.*exp(-landa_Ba133m*(deltat*j-deltat*1+deltat*c));


            end
            if (sum(temp_Xem127(i,:))<0.1*sum(temp_Xem127(d,:)) && f<1)
                f=i
                sum(temp_Xem127(i,:))
                0.1*sum(temp_Xem127(d,:))
            end
            c=c+1;
  %          sum(temp_Ba133m(i,:))
    end
           
            
   % else;
           % temp_C10t(i,:)=temp_C10t(i,:)+Y_C12_C10t.*landa_C10.*exp(-landa_C10*(t-1));
            
            
            
     
    temp_total2(i)=sum(sum(temp_Ba133(i,:))+sum(temp_Xe127(i,:))+sum(temp_Xe125(i,:))+sum(temp_Xe123(i,:))+sum(temp_Xe122(i,:))+sum(temp_Xem127(i,:)));
    temp_parcXem127(i)=sum(temp_Xem127(i,:));
    temp_parcXe127(i)=sum(temp_Xe127(i,:));
    temp_parcXe125(i)=sum(temp_Xe125(i,:));
    temp_parcXe123(i)=sum(temp_Xe123(i,:));
    temp_parcXe122(i)=sum(temp_Xe122(i,:));
    temp_parcBa133(i)=sum(temp_Ba133(i,:));
    temp_parcBa133m(i)=sum(temp_Ba133m(i,:));


    
    
end
    temp_total=temp_Xe127+temp_Xe125+temp_Xe123+temp_Xe122+temp_Xem127+temp_Ba133+temp_Ba133m;
    




    %% Dibuja la Actividad total en función del tiempo.
    figure;
    set(gca, 'FontSize', 16); 
    T=(0:t);
    hold on;
%    tottot=0.003122*0.38*temp_parcXem127+0.00285*0.68*temp_parcXem127+0.002535*0.1769*temp_parcBa133m;
%    plot(T*deltat,tottot/1000,'k','linewidth',2);
%    plot(T*deltat,(0.003122*0.38*temp_parcXem127)/1000,'r','linewidth',2);
%    plot(T*deltat,(0.00285*0.68*temp_parcXem127)/1000,'g','linewidth',2);
    plot(T*deltat,(temp_parcXem127)/1000,'r');
    plot(T*deltat,(temp_parcXem127)/1000,'g');
    plot(T*deltat,(temp_parcBa133m)/1000,'b')
    plot(T*deltat,(temp_parcXe127)/1000,'k')
    grid on
    title('Actividad/nC ');
    xlabel('Tiempo (s)');
    ylabel('Actividad (kBq)');
    legend('Total','Xe127*-172 keV','Xe127*-124 keV','Ba133*-275 keV','127xe');
 %   legend('Xe127*-172 keV','Xe127*-124 keV','Xe127','Xe125','Xe125','Xe123','Xe122');
%    axis([0 t*deltat 0 (0.68*max(temp_total2(:,1)/1000)+0.2*max(temp_total2(:,1)/1000))]); 
%     
max(temp_parcXem127)/1000
max(temp_parcXe127)/1000
max(temp_parcXe125)/1000
max(temp_parcXe123)/1000
max(temp_parcXe122)/1000
%% Actvidad 133Ba
    figure
    T=linspace(1,1.2614e9,40);
    act_Ba133t2=zeros(1,40);
    for i=1:(length(T))


    act_Ba133t2(i) =  landa_Ba133 .* sum(Y_Ba133) * pps * d .* exp(- landa_Ba133 * T(i));
    end
    
    
    TT=T/(24*365*3600);
    plot(TT,0.62*act_Ba133t2,'r','linewidth',2)
    title('Actividad residual CsI')
    xlabel('Tiempo (años)')
    ylabel('Actividad (Bq)')
    legend('133Ba-356 keV')

%% Actvidad 133mBa
    figure
    T=linspace(1,432000,120);
    act_Ba133mt2=zeros(1,120);
    for i=1:(length(T))

    act_Ba133mt2(i) =  landa_Ba133m .* sum(Y_Ba133m) * pps * d .* exp(- landa_Ba133m * T(i));
    
    end
    
    
    TT=T/(3600);
    plot(TT,0.1769*act_Ba133mt2,'r','linewidth',2)
    title('Actividad residual CsI')
    xlabel('Tiempo (horas)')
    ylabel('Actividad (Bq)')
    legend('133*Ba-275.925 keV')



%% Actvidad 133mBa
    figure
    T=linspace(1,10368000,120);
    act_Xe127t2=zeros(1,120);
    for i=1:(length(T))

    act_Xe127t2(i) =  landa_Xe127 .* sum(Y_Xe127) * pps * d .* exp(- landa_Xe127 * T(i));
    
    end
    
    
    TT=T/(24*3600);
    plot(TT,0.62*act_Xe127t2,'r','linewidth',2)
    title('Actividad residual CsI')
    xlabel('Tiempo (dias)')
    ylabel('Actividad (Bq)')
    legend('127Xe-356 keV')
% 
 %% Actividad total en función de tiempo
close all

