clear all; clc; close all;

%% Char. statyczna
clear all
h = 0.01;
stopTime = 20;
uInit = [-2:0.01:2];

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
% Inercja n-tego rzędu


