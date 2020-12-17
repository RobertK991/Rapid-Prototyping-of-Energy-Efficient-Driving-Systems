clear all; clc; close all;
load('K.mat');
load('nastawy.mat');
h = 0.01;
stopTime = 80;

model.licznik = 0.314656258763423;
model.mianownik = [1,1.00277037976700,0];

%% Char. statyczna
clear all
h = 0.01;
stopTime = 300;
uInit = [-2:0.1:2];  
% uInit = [-2:0.01:2];
for i=1:1:length(uInit)
    u_0 = uInit(i);
    sim('model_silnik_dykretny_lock.slx');
%     dane = y.data;
%     data = pi*(dane/180);       %Zmiana na rad - dodano do simulink
    yMean(i) = mean(y.data(end-100:end));
end

figure;
plot(uInit,yMean);
title('Static characteristic');
ylabel('y [rad]');
xlabel('u [V]');
grid on;

%% Punkt 3
h = 0.01;
u_0 = 1;
stopTime = 20;

sim('model_silnik_dykretny_lock.slx');
% wyjscie = pi*(y/180);       %Zmiana na rad - dodano do simulink
figure;
plot(y);
title(['Output u = ', num2str(u_0)]);
ylabel('y [rad]');
xlabel('u [V]');
grid on;
% Inercja pierwszego rzędu z całką

%% Identyfikacja - punkt 4
h = 0.01;
u_0 = 1;  % 50% 

w = numel(y.data);

data = y.data;
iden = ones(w,1)*u_0;

systemIdentification

%% K
% K = tf(tf2.Numerator, tf2.Denominator)

% Fit to estimation data: 99%, FPE: 0.000320976

% K = 
%            0.3147
%  -------------------------
%   s^2 + 1.003 s + 1.017e-12
%  
% Continuous-time transfer function.

% Pracujemy na tym:
% K =
%  
%      0.3147
%   -------------
%   s^2 + 1.003 s

%% Compare
figure;
plot(y,'r');
hold on;
step(K,20);

title('Compare')
legend('Plant','Model');
grid on;
ylabel('y [rad]');
xlabel('u [V]');
%% Punkt 5 i 6
% Regulator PID 
% Regulator P 
%Załadowanie danych z identyfikacji
load('K.mat');

h = 0.01;
uMax = 2;
iteracje = 4;
stopTime = 50;

for i = 0:iteracje
    u = (0.3 + 0.1*i)*uMax;
    [yStep, tStep] = step(K*u,0:h:stopTime);
    yMean = mean(yStep(end-100:end));
    eps = 0.5;
    
    %Metoda dwóch punktów
    indexMin = find(abs(yStep - 0.283*yMean)<eps);
    indexMin = min(indexMin);
    indexMax = find(abs(yStep - 0.632*yMean)<eps);
    indexMax = max(indexMax);
    
    % t1 i t2
    t1 = tStep(indexMax);
    t2 = tStep(indexMin);
    % Stała czasowa inercji I rzędu 
    T = 1.5*(t1 - t2);
    % Czas opóźnienia
    T0 = t1 - T;
    % Wzmocnienie modelu
    k = yMean/u;
    
    % Nastawy regulatorów QRD
    regP.kr(i+1) = T/(k*T0);
    regPID.kr(i+1) = 1.2 * (T/(k*T0));
    regPID.Ti(i+1) = 2*T0;
    regPID.Td(i+1) = 0.5*T0;
end

% Średnie nastawy regulatorów
regP.kr = mean(regP.kr);
regPID.kr = mean(regPID.kr);
regPID.Ti = mean(regPID.Ti);
regPID.Td = mean(regPID.Td);

% Dane modelu
model.licznik = 0.314656258763423;
model.mianownik = [1,1.00277037976700,0];
disp(['Nastawy regulatora P: kr = ', num2str(regP.kr)]);
disp(['Nastawy regulatora PID: kr = ', num2str(regPID.kr),' Ti = ',num2str(regPID.Ti),' Td = ',num2str(regPID.Td)]);

save('nastawy', 'regP', 'regPID')
%% Symulacja regulator P i PID

% Wartość zadana
u = 1;
stopTime = 80;

%Nastawy obliczone
load('nastawy.mat');

sub.F = 1;
sub.T = 1;
strojenieReczne('model7.slx', regP, regPID, sub,[],[]);


%% Punkt 8
close all;
stopTime = 80;
uInit = [-1, 1, pi];
set = [];   % puste - żeby wyśweitlić setValue
sub.F = 1;
sub.T = 1;
u = uInit(1);
kolory = ['b','r'];
legenda = [string(['Plant u = ',num2str(u)]),string(['Model u = ',num2str(u)])];
strojenieReczne('model7.slx', regP, regPID, sub, kolory, set, legenda);

set = 1;
u = uInit(2);
kolory = ['y','c'];
legenda = [string(['Plant u = ',num2str(u)]),string(['Model u = ',num2str(u)])];
strojenieReczne('model7.slx', regP, regPID, sub, kolory, set, legenda);

u = uInit(3);
kolory = ['m','g'];
legenda = [string(['Plant u = ',num2str(u)]),string(['Model u = ',num2str(u)])];
strojenieReczne('model7.slx', regP, regPID, sub, kolory, set, legenda);
%{
for i = 1:1:length(uInit)
u = uInit(i);
sub.F = 1;
sub.T = 1;
kolory = ['b','r'];

legenda = [string(['u = ',num2str(u)]),string(['u = ',num2str(u)])];
strojenieReczne('model7.slx', regP, regPID, sub, kolory, set, legenda);
set = 1;

% out = sim('model7.slx'); 
% 
% % Regulator P
% figure(1);
% subplot(length(uInit),1,i);
% yline(u,'--');
% hold on;
% plot(out.tout,out.yPlantP);
% plot(out.tout,out.yModelP);
% title(['Output, u = ', num2str(u), ', regulator P']);
% ylabel('y [rad]');
% xlabel('time');
% grid on;
% legend('Set value','Plant', 'Model');
% 
% % Regulator PID
% figure(2);
% subplot(length(uInit),1,i);
% yline(u,'--');
% hold on;
% plot(out.tout,out.yPlantPID);
% plot(out.tout,out.yModelPID);
% title(['Output, u = ', num2str(u), ', regulator PID']);
% ylabel('y [rad]');
% xlabel('time');
% grid on;
% legend('Set value','Plant', 'Model');

end
%}
%% Punkt 9 - dokończ go, - dwa optymalne regulatory - parametry
close all;
% krP = 0.954623487930613
% kr = 1.14554818551674
% ti = 6.267499999999997
% td = 1.566874999999999

% Ręczna zmiana nastaw 
stopTime = 90;
% wartośc zadana
u = 1;
%Nastawy reguatorów
krP = [1.2, 2];
krPID = [2, 4];
Ti = [ 7, 7];
Td = [0.1, 0.9];
%Zmiana częstotliwośći
freq = 0.1; 
% Do bloku chirp
freqInit = 0.05;
freqTime = 90;
freqTarget = 0.1;

licznik = length(krP);
for i = 1:1:licznik
 
% modelNazwa = 'model9_prost.slx';  %Wymuszenie prostokątnym
% modelNazwa = 'model9_sin.slx';    %Wymuszenie sinusem
modelNazwa = 'model9_sinZmienny.slx'; %Wym. sin z zmienną częstotliwością
% modelNazwa = 'model7.slx';        % Wymuszenie step

regP.kr = krP(i);
regPID.kr = krPID(i);
regPID.Ti = Ti(i);
regPID.Td = Td(i); 

sub.F = licznik;
sub.T = i;

out = strojenieReczne(modelNazwa, regP, regPID, sub);

% % Char bodego
% dataP = iddata(out.yPlantP,out.setValue,h);
% dataPID = iddata(out.yPlantPID,out.setValue,h);
% sysP = tfest(dataP,2)
% sysPID = tfest(dataPID,4)
% figure;
% bode(sysP);
% title(['kr = ', num2str(regP.kr), ', regulator P']);
% figure;
% bode(sysPID);
% title(['kr=', num2str(regPID.kr), ' Ti=', num2str(regPID.Ti), ' Td=', num2str(regPID.Td) ', regulator PID']);

end

%% Zmiany PID
%{
% close all;

u = 1;
freq = 0.1

krPID = [2];
% Ti = [1, 4, 6.26, 8];
Ti = [ 7];
Td = [0.1];

licznik = length(krPID);
licznikTi = length(Ti);
licznikTd = length(Td);
for i = 1:1:licznik
    for j = 1:1:licznikTi
        for k = 1:1:licznikTd
        
        regPID.kr = krPID(i);
        regPID.Ti = Ti(j);
        regPID.Td = Td(k); 

        out = sim('model9_prost.slx'); 
        out = sim('model9_prost.slx'); 
        
        % Regulator PID
        figure(i);
        subplot(licznik,1,k);   %zmiana j lub k 
        plot(out.tout, out.setValue,'k--');
        hold on;
        plot(out.tout,out.yPlantPID, 'b');
        plot(out.tout,out.yModelPID, 'r');
        title(['kr=', num2str(regPID.kr), ' Ti=', num2str(regPID.Ti), ' Td=', num2str(regPID.Td) ', regulator PID']);
        ylabel('y [rad]');
        xlabel('time');
        grid on;
        legend('Set value','Plant', 'Model');

        end
    end
end
%}

regP.kr = 1.2;  %Prostokątne: Przy większych wzmocnieniach są przeregulowania
%Sin: Im większe tym lepiej podąża za sinusem.
% regP.kr = 2;  % Większe przeregulowania, lepeij podąża za sinusem

%% Char bodego
% wy_1 = out(1);
% wy_2 = out(2);
% save('charBodeDane', 'wy_1', 'wy_2')
% load('charBodeDane.mat');

% Ręczna zmiana nastaw 
stopTime = 2000;
% wartośc zadana
u = 1;
% Do bloku chirp
freqInit = 0.05;
freqTime = 600;
freqTarget = 0.2;

regP.kr = (1.2 + 2)/2;

out = sim('model9_sinZmienny.slx'); 

figure(1);
plot(out.tout, out.setValue,'k--');
hold on;
plot(out.tout,out.yPlantP, 'b');
plot(out.tout,out.yModelP, 'r');
title(['kr = ', num2str(regP.kr), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');


%% 

% windowSize = 100; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;

[b,a] = butter(4,0.02,'low');
y = filter(b,a,out.yPlantP);

figure;
plot(y);
hold on;

% Znajduje piki
[pks, locs] = findpeaks(y);
%wykreśla piki

% scatter(locs, pks,'r');
% Piki > 0 
locs = locs(pks>0);
pks = pks(pks>0);
% figure;
scatter(locs, pks,'r');

bodeWar = 20*log10(pks);
% Bode z danych
figure;
plot(locs(1:end), bodeWar(1:end));
%bODE Z k
figure;
bode(K);

% Obie osie mają być w logorytmicznie !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%
dataP = iddata(out.yPlantP,out.setValue,h);
sysP = tfest(dataP,2)
figure;
bode(sysP);

%% Punkt 10 - sisotool

% sisotool(K)
load('sisoPID.mat');
%Sisotool - P
stopTime = 90;
u = 1;
regP.kr = 1.1437;

C = pid(C);
regPID.kr = C.Kp;
regPID.Ti = C.Ki;
regPID.Td = C.Kd; 

sub.F = 1;
sub.T = 1;
strojenieReczne('model7.slx', regP, regPID, sub);


%% 11 - dokończ robercik to !!!!!!!!

stopTime = 90;
u = 1;
% Siso tool
load('sisoNastawy.mat');
out = sim('model7.slx'); 

% Regulator P - sisotool
figure(1);
plot(out.tout, out.setValue,'k--');
hold on;
plot(out.tout,out.yPlantP);
plot(out.tout,out.yModelP);

% Regulator PID - sisotool
figure(2);
plot(out.tout, out.setValue,'k--');
hold on;
plot(out.tout,out.yPlantPID);
plot(out.tout,out.yModelPID);




% Regulator P Nastawy analityczne
load('nastawy.mat');
out = sim('model7.slx'); 
plot(out.tout,out.yPlantP);
plot(out.tout,out.yModelP);


% Title i opisy
title(['kr = ', num2str(regP.kr), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;

% Regulator PID
figure(2);
plot(out.tout, out.setValue,'k--');
hold on;
plot(out.tout,out.yPlantPID);
plot(out.tout,out.yModelPID);



% Title i opisy
title(['kr=', num2str(regPID.kr), ' Ti=', num2str(regPID.Ti), ' Td=', num2str(regPID.Td) ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;

figure(1);
legend('Set value','Plant', 'Model','Set value', 'SisotoolPlant', 'SisotoolModel');
figure(2);
legend('Set value','Plant', 'Model','Set value', 'SisotoolPlant', 'SisotoolModel');








