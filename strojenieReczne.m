function strojenieReczne(model, regP, regPID, sub)

out = sim(model); 

% Regulator P
figure(1);
subplot(sub.F,1,sub.T);
plot(out.tout, out.setValue,'k--');
hold on;
plot(out.tout,out.yPlantP, 'b');
plot(out.tout,out.yModelP, 'r');
title(['kr = ', num2str(regP.kr), ', regulator P']);
ylabel('y [rad]');
xlabel('time');
grid on;
legend('Set value','Plant', 'Model');


% Regulator PID
figure(2);
subplot(sub.F,1,sub.T);   %zmiana j lub k 
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