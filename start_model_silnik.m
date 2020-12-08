clear all; clc; close all;

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
K = tf(tf2.Numerator, tf2.Denominator)

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
step(tf1,20);

title('Compare')
legend('Plant','Model');
grid on;

%% Punkt 5 i 6
% Regulator PID 
% Regulator P 
%Załadowanie danych z identyfikacji
load('K.mat');

h = 0.01;
uMax = 2;
iteracje = 3;
stopTime = 50;

for i = 0:iteracje
    u = (0.3 + 0.2*i)*uMax;
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

%% Symulacja regulator P i PID

% Wartość zadana
u = 2.5;
stopTime = 80;

% Model i silnik - regulator P
out = sim('model7.slx'); 

%Wyświetlanie P
figure;
yline(u,'--');
hold on;
plot(out.tout,out.yPlantP);
plot(out.tout,out.yModelP);
title(['Output, u = ', num2str(u), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

% Regulator PID
% Wartość zadana

%Wyświetlanie PID
figure;
yline(u,'--');
hold on;
plot(out.tout,out.yPlantPID);
plot(out.tout,out.yModelPID);
title(['Output, u = ', num2str(u), ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

%% Punkt 8
stopTime = 80;
uInit = [-2, 2, pi];

for i = 1:1:length(uInit)
u = uInit(i);

out = sim('model7.slx'); 

% Regulator P
figure(1);
subplot(length(uInit),1,i);
yline(u,'--');
hold on;
plot(out.tout,out.yPlantP);
plot(out.tout,out.yModelP);
title(['Output, u = ', num2str(u), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

% Regulator PID
figure(2);
subplot(length(uInit),1,i);
yline(u,'--');
hold on;
plot(out.tout,out.yPlantPID);
plot(out.tout,out.yModelPID);
title(['Output, u = ', num2str(u), ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

end

%% Punkt 9 - dokończ go Robercik, - dwa optymalne regulatory - parametry
stopTime = 80;
% wartośc zadana
u = 2;


% Ręczna zmiana nastaw 
kr = [];
Ti = [];
Td = [];
licznik = length(kr);

for i = 1:1:licznik

out = sim('model7.slx'); 

% Regulator P
figure(1);
subplot(licznik,1,i);
yline(u,'--');
hold on;
plot(out.tout,out.yPlantP);
plot(out.tout,out.yModelP);
title(['Output, u = ', num2str(u), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

% Regulator PID
figure(2);
subplot(licznik,1,i);
yline(u,'--');
hold on;
plot(out.tout,out.yPlantPID);
plot(out.tout,out.yModelPID);
title(['Output, u = ', num2str(u), ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');

end

























