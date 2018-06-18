% C?digo de prueba de c?lculo de actividad con Zn
% (C) Daniel S?nchez Parcerisa 2018
% dsparcerisa@ucm.es

%% Par?metros a modificar:
clear all
dx = 0.1; % Espaciado de la malla (cm)
E0 = 140; % Energ?a inicial del haz (MeV)

%% Cargar datos
% Secciones de EXFOR.
load('CrossSections.mat');
C12_C10_CS=C12_C10_CS./1000;
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
landa_Bi206= log(2)/T_Bi206;
landa_Bi205= log(2)/T_Bi205;
landa_Bi204= log(2)/T_Bi204;
landa_Bi202= log(2)/T_Bi202;




% Natural abundances
Pb204_ab=0.014;
Pb206_ab=0.241;
Pb207_ab=0.221;
Pb208_ab=0.524;



% Cargar stopping powers (only for water, tissue, Zn, bone)
load('stoppingpowers.mat');

%% Fit secciones eficaces 206Pb_Bi
figure
Eval = 0:0.1:300; % MeV
hold off
hold on

plot(Pb206_202Bi_E,Pb206_202Bi_CS,'bo')
Pb206_202Bi_F = fit(Pb206_202Bi_E,Pb206_202Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb206_202Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb206_202Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb206_202Bi_F(Eval),'b-');

plot(Pb206_204Bi_E,Pb206_204Bi_CS,'go')
Pb206_204Bi_F = fit(Pb206_204Bi_E,Pb206_204Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb206_204Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb206_204Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb206_204Bi_F(Eval),'g-');

plot(Pb206_205Bi_E,Pb206_205Bi_CS,'ro')
Pb206_205Bi_F = fit(Pb206_205Bi_E,Pb206_205Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb206_205Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb206_205Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb206_205Bi_F(Eval),'r-');

axis([0 100 0 1.2]);


xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
legend('Bi202 data','Bi202 fit','Bi204 data','Bi204 fit','Bi205 data','Bi205 fit')
title('204 Pb cross sections')

%% Fit secciones eficaces 207Pb_Bi
figure
hold off
hold on

plot(Pb207_204Bi_E,Pb207_204Bi_CS,'bo')
Pb207_204Bi_F = fit(Pb207_204Bi_E,Pb207_204Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb207_204Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb207_204Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb207_204Bi_F(Eval),'b-');

plot(Pb207_205Bi_E,Pb207_205Bi_CS,'go')
Pb207_205Bi_F = fit(Pb207_205Bi_E,Pb207_205Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb207_205Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb207_205Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb207_205Bi_F(Eval),'g-');

plot(Pb207_206Bi_E,Pb207_206Bi_CS,'ro')
Pb207_206Bi_F = fit(Pb207_206Bi_E,Pb207_206Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb207_206Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb207_206Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb207_206Bi_F(Eval),'r-');




axis([0 100 0 1.2]);


xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
legend('Bi204 data','Bi204 fit','Bi205 data','Bi205 fit','Bi206 data','Bi206 fit')
title('206 Pb cross sections')



%% Fit secciones eficaces 208Pb_Bi
figure
hold off
hold on


plot(Pb208_205Bi_E,Pb208_205Bi_CS,'go')
Pb208_205Bi_F = fit(Pb208_205Bi_E,Pb208_205Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb208_205Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb208_205Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb208_205Bi_F(Eval),'g-');

plot(Pb208_206Bi_E,Pb208_206Bi_CS,'ro')
Pb208_206Bi_F = fit(Pb208_206Bi_E,Pb208_206Bi_CS,'smoothingspline','SmoothingParam',0.1);
Pb208_206Bi_F.p.coefs(1,:) = [0 0 0 0];
Pb208_206Bi_F.p.coefs(end,:) = [0 0 0 0];
plot(Eval,Pb208_206Bi_F(Eval),'r-');




axis([0 100 0 1.2]);


xlabel('Proton energy (MeV)')
ylabel('Cross section (barn)')
legend('Bi205 data','Bi205 fit','Bi206 data','Bi206 fit')





















%% Create fit for stopping power
figure


S_pb_F = fit(E_keV_PMMA,S_Pb,'smoothingspline','SmoothingParam',0.002)
S_pb_F.p.coefs(1,:) = [0 0 0 0];
S_pb_F.p.coefs(end,:) = [0 0 0 0];

loglog(E_keV_PMMA,S_Pb,'go')
hold on;
loglog(E_keV_PMMA,S_pb_F(E_keV_PMMA),'g-')
legend('Pb data', 'Fit Pb');
xlabel('Proton energy (MeV)');
ylabel('Stopping power (MeV/(cm2/mg))');
%%
close all

%% Calcular (sin straggling)

AvNmbr = 6.022140857e23;
rho_pb=11.34; %g/cm3
rho_pb_mol=207.19;%g/moñ
rho_pb_A=AvNmbr*rho_pb/rho_pb_mol; %at/cm3
pps=1; %proton per second

x = 0:dx:4; % posiciones en cm.
Epb = nan(size(x));

Ddeppb = zeros(size(Epb));

currentEpb = E0;




%PLOMO

Y_Pb206_Bi202=zeros(size(x));
Y_Pb206_Bi204=zeros(size(x));
Y_Pb206_Bi205=zeros(size(x));
Y_Pb207_Bi204=zeros(size(x));
Y_Pb207_Bi205=zeros(size(x));
Y_Pb207_Bi206=zeros(size(x));
Y_Pb208_Bi205=zeros(size(x));
Y_Pb208_Bi206=zeros(size(x));



for i=1:(numel(x)-1)

    S_pb = max(0,1000*S_pb_F(currentEpb*1000)); % MeV/(g/cm2)
    

    S1pb = (S_pb*rho_pb); % MeV/cm
    
 
    deltaEpb = dx*S1pb; %MeV
    


    currentEpb = currentEpb - deltaEpb; %MeV
    

    S_pb2 = max(0,1000*S_pb_F(currentEpb*1000));% MeV/(g/cm2)
    

    S2pb = (S_pb2*rho_pb); % MeV/cm3

    Ddeppb(i) = deltaEpb; %MeV
    
    % Zn part

    E1pb = Epb(i);

    E2pb = currentEpb;

    

    
    % Tissue (simplified only)
    sigma_206_202 = 0.5 * (max(0,Pb206_202Bi_F(E1pb)) + max(0,Pb206_202Bi_F(E2pb)));
    sigma_206_204 = 0.5 * (max(0,Pb206_204Bi_F(E1pb)) + max(0,Pb206_204Bi_F(E2pb)));
    sigma_206_205 = 0.5 * (max(0,Pb206_204Bi_F(E1pb)) + max(0,Pb206_205Bi_F(E2pb)));
    sigma_207_204 = 0.5 * (max(0,Pb207_204Bi_F(E1pb)) + max(0,Pb207_204Bi_F(E2pb)));
    sigma_207_205 = 0.5 * (max(0,Pb207_205Bi_F(E1pb)) + max(0,Pb207_205Bi_F(E2pb)));
    sigma_207_206 = 0.5 * (max(0,Pb207_206Bi_F(E1pb)) + max(0,Pb207_206Bi_F(E2pb)));
    sigma_208_205 = 0.5 * (max(0,Pb208_205Bi_F(E1pb)) + max(0,Pb208_205Bi_F(E2pb)));
    sigma_208_206 = 0.5 * (max(0,Pb208_206Bi_F(E1pb)) + max(0,Pb208_206Bi_F(E2pb)));
    
    Y_Pb206_Bi202(i) =pps * rho_pb_A *Pb206_ab * sigma_206_202 * 1e-24 * dx;
    Y_Pb206_Bi204(i)=pps * rho_pb_A *Pb206_ab * sigma_206_204 * 1e-24 * dx;
    Y_Pb206_Bi205(i)=pps * rho_pb_A *Pb206_ab * sigma_206_205 * 1e-24 * dx;
    Y_Pb207_Bi204(i)=pps * rho_pb_A *Pb207_ab * sigma_207_204 * 1e-24 * dx;
    Y_Pb207_Bi205(i)=pps * rho_pb_A *Pb207_ab * sigma_207_205 * 1e-24 * dx;
    Y_Pb207_Bi206(i)=pps * rho_pb_A *Pb207_ab * sigma_207_206 * 1e-24 * dx;
    Y_Pb208_Bi205(i)=pps * rho_pb_A *Pb207_ab * sigma_208_205 * 1e-24 * dx;
    Y_Pb208_Bi206(i)=pps * rho_pb_A *Pb207_ab * sigma_208_205 * 1e-24 * dx;

    
    
end

%% Create plots

% Figure in PMMA
figure
%subplot(2,1,1)
%plot(x,Eb)
yyaxis right
xlabel('Depth (cm)');
ylabel('Dose (a.u.)')
title('Plomo');
hold on
plot(x,100*Ddeppb)
legend('Dose')
set(gca,'FontSize',14)
axis([0 10 0 (max(100*Ddeppb)+50)]);
%subplot(2,1,2)
yyaxis left
grid on
%title('Yields of different species (per incoming proton)');
hold on
Y_202pb = Y_Pb206_Bi202;
Y_204pb = Y_Pb206_Bi204+Y_Pb207_Bi204;
Y_205pb = Y_Pb206_Bi205+Y_Pb207_Bi205+Y_Pb208_Bi205;
Y_206pb = Y_Pb207_Bi205+Y_Pb208_Bi205;
plot(x,Y_202pb,'k'); hold on
plot(x,Y_204pb,'c')
plot(x,Y_205pb,'m')
plot(x,Y_206pb,'y');
legend('Bi202','Bi204','Bi204','Bi206','Location', 'northwest');
axis([0 3 0 max(max(Y_205pb,Y_206pb)+0.2*max(Y_205pb,Y_206pb))]);
xlabel('Depth (cm)');
ylabel('Yield');
set(gca,'FontSize',14)

%% PET activity (t)
calcTimes = [0 3600 10000 20000]; % s


act_202pb = zeros(numel(calcTimes), numel(x));
act_204pb = zeros(numel(calcTimes), numel(x));
act_205pb = zeros(numel(calcTimes), numel(x));
act_206pb = zeros(numel(calcTimes), numel(x));

for i=1:numel(calcTimes)
    
     % Bone
    act_202pb(i,:) = landa_Bi202 .* Y_202pb .* exp(- landa_Bi202 * calcTimes(i));
    act_204pb(i,:) = landa_Bi204 .* Y_204pb .* exp(- landa_Bi204 * calcTimes(i));    
    act_205pb(i,:) = landa_Bi205 .* Y_205pb .* exp(- landa_Bi205 * calcTimes(i));
    act_206pb(i,:) = landa_Bi206 .* Y_206pb .* exp(- landa_Bi206 * calcTimes(i)); 
    
end
act_total = (0.154*act_202pb+0.0051*act_204pb+0.0041*act_205pb+0.00003*act_206pb);

G=figure;
I=figure;
AAA=length(calcTimes);

for i=1:numel(calcTimes)
 
    figure(I);
    subplot(2,2,i)
    yyaxis left
    plot(x, act_total(i,:),'b')
    hold on
    plot(x, 0.154*act_202pb(i,:),'r')
    plot(x, 0.0051*act_204pb(i,:),'y')
    plot(x, 0.0041*act_205pb(i,:),'c')
    plot(x, 0.00003*act_206pb(i,:),'g')
    title(sprintf('Activity at t=%i s (511 keV) ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Beta+ emitters / s ');
    legend('Total','Bi202','Bi204','Bi205','Bi206', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,100*Ddeppb);
    ylabel('-dE/dx')
    axis([0 3 0 max(100*Ddeppb)+0.2*max(100*Ddeppb)]);
end

for i=1:numel(calcTimes)
 
    figure(G);
    subplot(2,2,i)
    yyaxis left
    hold on
    plot(x, 0.82*act_204pb(i,:),'r')
    plot(x, 0.84*act_202pb(i,:),'y')
    plot(x, 0.408*act_206pb(i,:),'c')
    plot(x, 0.305*act_206pb(i,:),'m')
    plot(x, 0.61*act_202pb(i,:),'g')
    plot(x, 0.311*act_205pb(i,:),'b')
    plot(x, 0.99*act_206pb(i,:),'k')
    plot(x, 0.66*act_206pb(i,:),'r--')
    plot(x, 0.99*act_204pb(i,:),'k--')
    plot(x, 0.99*act_202pb(i,:),'m--')
    plot(x, 0.59*act_204pb(i,:),'b--')
    plot(x, 0.325*act_205pb(i,:),'k--')
    plot(x, 0.32*act_206pb(i,:),'b-.')
    

    title(sprintf('Activity at t=%i s (gammas)  ',calcTimes(i)));
    xlabel('Depth (cm)')
    ylabel('Emision / s ');
    legend('375 keV','422','516','537','658','703','803','881','900','960','984','1764','1781', 'Location', 'northwest');
    set(gca, 'FontSize', 16)   
  
    
    yyaxis right
    grid on
    plot(x,100*Ddeppb);
    ylabel('-dE/dx')
    axis([0 3 0 max(100*Ddeppb)+0.2*max(100*Ddeppb)]);
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
a=600;
c=1;
t=3600;
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
    plot(T*deltat,(temp_totalp2(:,1)));
    hold on;
    plot(T*deltat,temp_parcO15p(:,1));
    plot(T*deltat,temp_parcC11p(:,1));
    plot(T*deltat,temp_parcN13p(:,1));
    plot(T*deltat,temp_parcC10p(:,1));
    title('Actividad total en a lo largo del tiempo PMMA');
    xlabel('Tiempo (s)');
    ylabel('Beta+ emitters / s ');
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


