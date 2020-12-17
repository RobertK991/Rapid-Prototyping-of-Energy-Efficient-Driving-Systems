function [out] = strojenieReczne(model, regP, regPID, sub, kolor, setVal, legenda)
% strojenieReczne(model, regP, regPID, sub, kolor, setVal)
% model     - Nazwa modelu do urochomienia
% reP       - nastawy regulatora P
% regPID    - nastawy regulatora PID
% sub       - sub.F - ile ma być obrazków w oknie, sub.T - który to ma być
%             obrazek
% kolor     - kolory do wykresów - tylko podstawowe matlaba w 'r' formacie
% setVal    - czy ma być wykreślana setValue
% legenda   - należy dodać legendę w "..." 

if (~exist('kolor','var') || isempty(kolor))
    kolor = ['b', 'r'];
end

if (~exist('legend') || isempty(legenda))
    legenda = ["Plant", "Model"];
end

out = sim(model); 
set_figure_toscreen(1)
% Regulator P
figure(1);
subplot(sub.F,1,sub.T);
if (~exist('setVal','var') || isempty(setVal))
    plot(out.tout, out.setValue,'k--','DisplayName','Set value');
end
hold on;
plot(out.tout,out.yPlantP, kolor(1),'DisplayName',legenda(1));
plot(out.tout,out.yModelP, kolor(2),'DisplayName',legenda(2));
% plot(out.tout,out.yPlantP,'r');
% plot(out.tout,out.yModelP,'b');
title(['kr = ', num2str(regP.kr), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend show;
% Regulator PID
figure(2);
subplot(sub.F,1,sub.T);   %zmiana j lub k 
if (~exist('setVal','var') || isempty(setVal))
    plot(out.tout, out.setValue,'k--','DisplayName','Set value');
end
hold on;
plot(out.tout,out.yPlantPID, kolor(1),'DisplayName',legenda(1));
plot(out.tout,out.yModelPID, kolor(2),'DisplayName',legenda(2));
% plot(out.tout,out.yPlantPID,'r');
% plot(out.tout,out.yModelPID,'b');
title(['kr=', num2str(regPID.kr), ' Ti=', num2str(regPID.Ti), ' Td=', num2str(regPID.Td) ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;

legend show;
end