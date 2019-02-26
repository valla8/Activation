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
dx=0.002;      %Paso del intervalo (cm)
xref=1.0;       %Distancia que va a simular, poner un número acorde a la energia inicial.
E0=9;        %Energía inicial del haz

pps=5e10; %protones/segundo
MeVJ=1.6e-13;
landa_F18 =  log(2) / 6586;
landa_Ca39 = log(2) /0.859;
O18_fraction=0.00;
%% Calcular (sin straggling)

%CsI=NaCl

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
rho_NaCl=2.16;
rho_KCl=1.987;

W_ele=[ 39.1 35.4 ];

Comp_KCl=[0.5 0.5]

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
rho_KCl_A =  rho_KCl * AvNmbr / sum(Comp_KCl.*W_ele); % molecules / cm3

%Calculamos la densidad de cada isótopo multiplicando por su peso y su
%abundancia. Se hace para cada material bone(b) tissue(t) PMMA(p)
%adipose(a) water ()
rho_K_A = rho_KCl_A * Comp_KCl(1); % atoms/cm3
rho_Cl_A = rho_KCl_A * Comp_KCl(2); % atoms/cm3

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
    S_KCl = max(0,1000*S_KCl_F(currentE*1000)); % MeV/(g/cm2)
    
    %Multiplicamospor la densidad al poder de frenado
    S1 = S_KCl*rho_KCl; % MeV/cm
    
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
    Y_39Ca = zeros(size(x));
    

    
     for i=1:length(x)

     
         
    
    % Water (full + simplified)

   

     Y_39Ca(i) = Y_39Ca(i)  +  rho_K_A * max(0,K39_Ca39_F(E(i))) * 1e-24 * dx;
    
%     EE(i) = EE(i) + Fl(k) * Ddep(k);
    
    
    
     end
        
     



    
    
    % end
%%




figure
hold on
plot(E,K39_Ca39_F(E))
legend('Ca39')




%%
figure
hold on

 plot(x,Y_39Ca)
legend('39Ca')

 %% PET activity (t) para un cierto número de protones que llegan simultaneamente
% %Valores para los que se pretende calcular la actividad
calcTimes = [0 30 60 120 360 ]; % s

%Definimos las variables
act_Ca39 = zeros(numel(calcTimes), numel(x));

nC=1/1.6e-10;




%Calculo Actividad
%Como hemos introducido el p/s en el Yield roduction tenemos que calcular
%la actividad del número de protones en el intervalo de tiempo que estamos simulando.
for i=1:numel(calcTimes)
    % Water
    act_Ca39(i,:) =  landa_Ca39 .* Y_39Ca .* exp(- landa_Ca39 * calcTimes(i));
    

   
end
act_total = (act_Ca39);
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
act_Ca39t = zeros(1,360);
TT=[1:360];

for i=1:360


    act_Ca39t(i) =  landa_Ca39 .* sum(Y_39Ca) .* exp(- landa_Ca39 * i);
    
end


plot(TT,act_Ca39t)
hold on
legend('39Ca')


%% Actividad con el tiempo
%Aquí se calcula la actividad en función del tiempo cuando el tiempo de
%irradiación es superior al delta de t, es decir, un caso clínico real
%usando un ciclotrón (no está preparado todavía para la simualción de
%pulsos de un sincrotrón). Los parámetros de tiempo de irradiación y de
%decaimiento se introducen al principio

deltat=0.01;      %Inervalo de tiempo de las simulaciones
a=10/deltat;  %Tiempo de irradación del haz (s)
t=30/deltat;  %Tiempo total de la simulación

c=1;
temp_Ca39=zeros(t+1,numel(x));

temp_total2=zeros(t+1,1);
temp_parcCa39=zeros(t+1,1);

%CALCULO
%Para hacer el calculo suponemos que cada intervalod de tiempo llegan un
%número de protones que se frenan inmediatamente y genera los isótopos
%correspondientes. En cada iteración se calcula el número de los isótopos
%de las iteraciones anteriores que se ha desintegrado y además si estamos
%en el tiempo de irradiación se suma la contribución de los protones que
%llegan.
d=0
buff=0
buff1=sum(Y_39Ca);
for i=1:(t+1)
    b=i;
    if buff<buff1
           d=d+1;
    %if i<a
        for j=1:b
            
            temp_Ca39(i,:)=temp_Ca39(i,:)+pps*deltat * Y_39Ca.*landa_Ca39.*exp(-landa_Ca39*(deltat*j-deltat*1));
            
        end
           buff=sum(temp_Ca39(i,:))/pps;
    else
            for j=1:a;
                
            temp_Ca39(i,:)=temp_Ca39(i,:)+pps*deltat * Y_39Ca.*landa_Ca39.*exp(-landa_Ca39*(deltat*j-deltat*1+deltat*c));


            
            end
            if (sum(temp_Ca39(i,:))<0.1*sum(temp_Ca39(d,:)) && f<1)
                f=i
                sum(temp_Ca39(i,:))
                0.1*sum(temp_Ca39(d,:))
            end
            c=c+1;
    end
           
            
   % else;
           % temp_C10t(i,:)=temp_C10t(i,:)+Y_C12_C10t.*landa_C10.*exp(-landa_C10*(t-1));
            
            
            
     
    temp_total2(i)=sum(sum(temp_Ca39(i,:)));
    temp_parcCa39(i)=sum(temp_Ca39(i,:));


    
    
end
    temp_total=temp_Ca39;
    




    %% Dibuja la Actividad total en función del tiempo.
    figure;
    set(gca, 'FontSize', 16); 
    T=(0:t);
    hold on;
    plot(T*deltat,0.005*(temp_parcCa39)/1000,'r');
    grid on
    title('Actividad total en a lo largo del tiempo');
    xlabel('Tiempo (s)');
    ylabel('Actividad (kBq)');
    legend('Ca39');
 %   axis([ 0 10 0 300 ])
 %   legend('Xe127*-172 keV','Xe127*-124 keV','Xe127','Xe125','Xe125','Xe123','Xe122');
%    axis([0 t*deltat 0 (0.68*max(temp_total2(:,1)/1000)+0.2*max(temp_total2(:,1)/1000))]); 
%     

% 
% %% Actividad total en función de tiempo
% %Igual que el apartado anterior pero con otra fórmula. Es una comprobación
% %de que estaba bien hecho.
% cont1=0;
% cont2=0;
% c=1;
% int_totalp2=zeros(t+1,1);
% int_totalt2=zeros(t+1,1);
% int_parcC11t=zeros(t+1,1);
% int_parcO15t=zeros(t+1,1);
% int_parcN13t=zeros(t+1,1);
% int_parcC10t=zeros(t+1,1);
% int_parcC11p=zeros(t+1,1);
% int_parcO15p=zeros(t+1,1);
% int_parcN13p=zeros(t+1,1);
% int_parcC10p=zeros(t+1,1);
% int_C11t=zeros(t+1,numel(x));
% int_O15t=zeros(t+1,numel(x));
% int_N13t=zeros(t+1,numel(x));
% int_C10t=zeros(t+1,numel(x));
% int_C11p=zeros(t+1,numel(x));
% int_O15p=zeros(t+1,numel(x));
% int_N13p=zeros(t+1,numel(x));
% int_C10p=zeros(t+1,numel(x));
% for i=2:t+1;
%     b=i;
%     if i<a
%         cont1=cont1+1;
%         for j=1:b;
%             
%             %Tissue
%             int_C10t(i,:)=int_C10t(i,:)+deltat*Y_C10t.*(-exp(-landa_C10*(deltat*j))-exp(-landa_C10*(deltat*j-deltat*1)));
%             int_C11t(i,:)=int_C11t(i,:)+deltat*Y_C11t.*(-exp(-landa_C11*(deltat*j))+exp(-landa_C11*(deltat*j-deltat*1)));
%             int_N13t(i,:)=int_N13t(i,:)+deltat*Y_N13t.*(-exp(-landa_N13*(deltat*j))+exp(-landa_N13*(deltat*j-deltat*1)));;
%             int_O15t(i,:)=int_O15t(i,:)+deltat*Y_O15t.*(-exp(-landa_O15*(deltat*j))+exp(-landa_O15*(deltat*j-deltat*1)));;
%             
%             %PMMA
%             int_C10p(i,:)=int_C10p(i,:)+deltat*Y_C10p.*(-exp(-landa_C10*(deltat*j))+exp(-landa_C10*(deltat*j-deltat*1)));
%             int_C11p(i,:)=int_C11p(i,:)+deltat*Y_C11p.*(-exp(-landa_C11*(deltat*j))+exp(-landa_C11*(deltat*j-deltat*1)));
%             int_N13p(i,:)=int_N13p(i,:)+deltat*Y_N13p.*(-exp(-landa_N13*(deltat*j))+exp(-landa_N13*(deltat*j-deltat*1)));
%             int_O15p(i,:)=int_O15p(i,:)+deltat*Y_O15p.*(-exp(-landa_O15*(deltat*j))+exp(-landa_O15*(deltat*j-deltat*1)));
%            
%         end
% 
%     else
%        
%            cont2=cont2+1;
%         for j=1:a;
%             
%             %Tissue
%             int_C10t(i,:)=int_C10t(i,:)+deltat*Y_C10t.*(-exp(-landa_C10*(deltat*j+deltat*c))+exp(-landa_C10*(deltat*j-deltat*1+deltat*c)));
%             int_C11t(i,:)=int_C11t(i,:)+deltat*Y_C11t.*(-exp(-landa_C11*(deltat*j+deltat*c))+exp(-landa_C11*(deltat*j-deltat*1+deltat*c)));
%             int_N13t(i,:)=int_N13t(i,:)+deltat*Y_N13t.*(-exp(-landa_N13*(deltat*j+deltat*c))+exp(-landa_N13*(deltat*j-deltat*1+deltat*c)));
%             int_O15t(i,:)=int_O15t(i,:)+deltat*Y_O15t.*(-exp(-landa_O15*(deltat*j+deltat*c))+exp(-landa_O15*(deltat*j-deltat*1+deltat*c)));
%             
%             %PMMA
%             int_C10p(i,:)=int_C10p(i,:)+deltat*Y_C10p.*(-exp(-landa_C10*(deltat*j+deltat*c))+exp(-landa_C10*(deltat*j-deltat*1+deltat*c)));
%             int_C11p(i,:)=int_C11p(i,:)+deltat*Y_C11p.*(-exp(-landa_C11*(deltat*j+deltat*c))+exp(-landa_C11*(deltat*j-deltat*1+deltat*c)));
%             int_N13p(i,:)=int_N13p(i,:)+deltat*Y_N13p.*(-exp(-landa_N13*(deltat*j+deltat*c))+exp(-landa_N13*(deltat*j-deltat*1+deltat*c)));
%             int_O15p(i,:)=int_O15p(i,:)+deltat*Y_O15p.*(-exp(-landa_O15*(deltat*j+deltat*c))+exp(-landa_O15*(deltat*j-deltat*1+deltat*c)));
%            
%         end
%             c=c+1;
%             
% 
%       
%     end
%            
%     int_totalt2(i)=sum(int_C10t(i,:))+sum(int_C11t(i,:))+sum(int_N13t(i,:))+sum(int_O15t(i,:));
%     int_parcC11t(i)=sum(int_C11t(i,:));
%     int_parcO15t(i)=sum(int_O15t(i,:));
%     int_parcN13t(i)=sum(int_N13t(i,:));
%     int_parcC10t(i)=sum(int_C10t(i,:));
%     
%     int_totalp2(i)=sum(int_C10p(i,:))+sum(int_C11p(i,:))+sum(int_N13p(i,:))+sum(int_O15p(i,:));
%     int_parcC11p(i)=sum(int_C11p(i,:));
%     int_parcO15p(i)=sum(int_O15p(i,:));
%     int_parcN13p(i)=sum(int_N13p(i,:));
%     int_parcC10p(i)=sum(int_C10p(i,:));
% 
%     
% end
% 
% %% Actividadd total 2
% 
%     figure;
%     plot(T*deltat,(int_totalp2(:,1)));
%     hold on
%     plot(T*deltat,int_parcO15p(:,1))
%     plot(T*deltat,int_parcC11p(:,1))
%     plot(T*deltat,int_parcC10p(:,1))
%     plot(T*deltat,int_parcN13p(:,1))
%     title('Actividad total en a lo largo del tiempo PMMA');
%     xlabel('Tiempo (s)')
%     ylabel('Beta+ emitters / s ');
%     legend('Actividad total','O15','C11','C10','N13','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     axis([0 t*deltat 0 (max(int_totalp2(:,1))+0.2*max(int_totalp2(:,1)))]);
% 
%     
%     
% %% Actividad durante Bean-On
% %Aquí calculamos la cantidad de pares de 511 keV que se producen durante la
% %irradiación y los que ocurren después en función del espacio, hasta el valor del tiempo en
% %segundos que hay de la variable de abajo
% 
% beam=zeros(1,length(x));
% beam_onp=zeros(1,length(x));
% beam_offp=zeros(1,length(x));
% C10_onp=zeros(1,length(x));
% C10_offp=zeros(1,length(x));
% C11_onp=zeros(1,length(x));
% C11_offp=zeros(1,length(x));
% N13_onp=zeros(1,length(x));
% N13_offp=zeros(1,length(x));
% O15_onp=zeros(1,length(x));
% O15_offp=zeros(1,length(x));
% 
% for i=1:tt
%     beam=beam+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
%     if i<a
%         beam_onp=beam_onp+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
%         C10_onp=C10_onp+int_C10p(i,:);
%         C11_onp=C11_onp+int_C11p(i,:);
%         N13_onp=N13_onp+int_N13p(i,:);
%         O15_onp=O15_onp+int_O15p(i,:);
%     else
%         beam_offp=beam_offp+int_C10p(i,:)+int_C11p(i,:)+int_N13p(i,:)+int_O15p(i,:);
%         C10_offp=C10_offp+int_C10p(i,:);
%         C11_offp=C11_offp+int_C11p(i,:);
%         N13_offp=N13_offp+int_N13p(i,:);
%         O15_offp=O15_offp+int_O15p(i,:);
%     end
% end
% 
% figure;
% 
% %subplot(2,2,1);
%     yyaxis right
%     grid on
%     plot(x,100*Ddepp,'k--');
%     ylabel('-dE/dx')
%     axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
%     yyaxis left
% hold on;
% plot(x,beam_onp(1,:),'b');
% plot(x,beam_offp(1,:),'r');
%     title('Actividad total ');
%     xlabel('z (cm)')
%     ylabel('Beta+ emitters  ');
%     legend('Beam On','Beam Off','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     grid on;
% [f,g]=min(Ddepp);
% axis([ 0 (ceil(g*dx)) 0 (max(beam_offp)+0.2*max(beam_offp))]);
% 
% figure;
% %subplot(2,2,2);
%     yyaxis right
%     grid on
%     plot(x,100*Ddepp,'k--');
%     ylabel('-dE/dx')
%     axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
%     yyaxis left
% hold on;
% plot(x,C11_onp(1,:),'b');
% plot(x,C11_offp(1,:),'r');
%     title('C11 ');
%     xlabel('z (cm)')
%     ylabel('Beta+ emitters ');
%     legend('Beam On','Beam Off','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     grid on;
% [f,g]=min(Ddepp);
% axis([ 0 (ceil(g*dx)) 0 (max(C11_offp)+0.2*max(C11_offp))]);
% 
% figure;
% %subplot(2,2,3);
%     yyaxis right
%     grid on
%     plot(x,100*Ddepp,'k--');
%     ylabel('-dE/dx')
%     axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
%     yyaxis left
% hold on;
% plot(x,O15_onp(1,:),'b');
% plot(x,O15_offp(1,:),'r');
%     title('O15 ');
%     xlabel('z (cm)')
%     ylabel('Beta+ emitters  ');
%     legend('Beam On','Beam Off','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     grid on;
% [f,g]=min(Ddepp);
% axis([ 0 (ceil(g*dx)) 0 (max(O15_offp)+0.2*max(O15_offp))]);
% 
% figure;
% %subplot(2,2,4);
%     yyaxis right
%     grid on
%     plot(x,100*Ddepp,'k--');
%     ylabel('-dE/dx')
%     axis([0 17 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
%     yyaxis left
% hold on;
% plot(x,N13_onp(1,:),'b');
% plot(x,N13_offp(1,:),'r');
%     title('N13 ');
%     xlabel('z (cm)')
%     ylabel('Beta+ emitters  ');
%     legend('Beam On','Beam Off','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     grid on;
% [f,g]=min(Ddepp);
% axis([ 0 (ceil(g*dx)) 0 (max(N13_offp)+0.2*max(N13_offp))]);
% 
% figure
%     yyaxis right
%     grid on
%     plot(x,100*Ddepp,'k--');
%     ylabel('-dE/dx')
%     axis([0 5 0 max(100*Ddepp)+0.3*max(100*Ddepp)]);
%     yyaxis left
% hold on;
% plot(x,beam(1,:),'b');
%     title('Total ');
%     xlabel('z (cm)')
%     ylabel('Beta+ emitters  ');
%     legend('Beam On','Dose','Location', 'northeast');
%     set(gca, 'FontSize', 16) 
%     grid on;
% [f,g]=min(Ddepp);
% axis([ 0 (ceil(g*dx)) 0 (max(beam(1,:))+0.2*max(beam(1,:)))]);
% 
% 
% 
% 
% 
