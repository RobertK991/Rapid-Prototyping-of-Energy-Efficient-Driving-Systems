function [out] = strojenieReczne(model, regP, regPID, sub, kolor, setVal, legenda, stylLini)
% [out] = strojenieReczne(model, regP, regPID, sub, kolor, setVal, legenda, stylLini)
% model     - Nazwa modelu do urochomienia
% reP       - nastawy regulatora P
% regPID    - nastawy regulatora PID
% sub       - sub.F - ile ma być obrazków w oknie, sub.T - który to ma być
%             obrazek
% kolor     - kolory do wykresów - w rgb [0 0 0]
% setVal    - czy ma być wykreślana setValue
% legenda   - należy dodać legendę w "..." 
% stylLini  - styl linii do wykresów, w "...", domyślnie "-"

% Sprawdz czy istnieje zmienna kolor, jeśli nie to nadaje kolor domyślny
if (~exist('kolor','var') || isempty(kolor))
    kolor = [0 0 1, 1 0 0];
end
% Sprawdz czy istnieje zmienna legenda, jeśli nie to nadaje legende
% domyślną
if (~exist('legenda') || isempty(legenda))
    legenda = ["Plant", "Model"];
end
% Sprawdz czy istnieje zmienna stylLini, jeśli nie to nadaje stylLini domyślny
if (~exist('stylLini') || isempty(stylLini))
    stylLini = ["-", "--"];
end
% Uruchamia model
out = sim(model); 
% Wybiera ekran na którym ma się wyświetlać.                            !!!
set_figure_toscreen(1)
% Regulator P
figure(1);
subplot(sub.F,1,sub.T);
% Sprawdza czy wyświetlić setValue
if (exist('setVal','var') && ~isempty(setVal))
    plot(out.tout, out.setValue,'k:','DisplayName','Set value');
end
hold on;
% Wyświetla dane z simulinka
plot(out.tout,out.yPlantP, 'Color',kolor(1:3),'DisplayName',legenda(1),'LineStyle', stylLini(1));
plot(out.tout,out.yModelP, 'Color',kolor(4:6),'DisplayName',legenda(2),'LineStyle', stylLini(2));
% Tytuły okien
title(['kr = ', num2str(regP.kr), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
% Pokazuje legendę
legend show;
% Regulator PID
figure(2);
subplot(sub.F,1,sub.T);   %zmiana j lub k 
% Sprawdza czy wyświetlić setValue
if (exist('setVal','var') && ~isempty(setVal))
    plot(out.tout, out.setValue,'k:','DisplayName','Set value');
end
hold on;
% Wyświetla dane z simulinka dla PID
plot(out.tout,out.yPlantPID, 'Color',kolor(1:3),'DisplayName',legenda(1),'LineStyle', stylLini(1));
plot(out.tout,out.yModelPID, 'Color',kolor(4:6),'DisplayName',legenda(2),'LineStyle', stylLini(2));
% Tytuły okien
title(['kr=', num2str(regPID.kr), ' Ti=', num2str(regPID.Ti), ' Td=', num2str(regPID.Td) ', regulator PID']);
ylabel('y [rad]');
xlabel('time');
grid on;
% Pokazuje legendę
legend show;
end