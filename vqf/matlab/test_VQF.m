close all; clear all; clc;

rawDataHuro = readtable("log\mag_biases\Huro_06_29_1533.csv");
% rawDataXsen = readtable("log\XSens_06_23_1413.csv");
rawGry = [rawDataHuro.gx, rawDataHuro.gy, rawDataHuro.gz];
rawAcc = [rawDataHuro.ax, rawDataHuro.ay, rawDataHuro.az];
rawMag = [rawDataHuro.mx, rawDataHuro.my, rawDataHuro.mz];

param.magNormTh = 0.1;

rad2deg = 180/pi;
acc_ts = 1/(1140/2);
gry_ts = 1/(1140/2);
mag_ts = 0.01;
filter_VQF = VQF(gry_ts, acc_ts, mag_ts);
output = filter_VQF.updateBatch(rawGry, rawAcc, rawMag);

time_huro = acc_ts : acc_ts : acc_ts*length(output.quat6D);
% time_xsen = 0.01 : 0.01 : 0.01*length(rawDataXsen.Time);

magNormDip = NaN([length(output.quat6D), 2]);
for idx = 1:length(output.quat6D)
    magEarth = quatRotate(output.quat6D(idx,:), [rawDataHuro.mx(idx), rawDataHuro.my(idx), rawDataHuro.mz(idx)]);
    magNormDip(idx,1) = norm(magEarth);
    magNormDip(idx,2) = -asin(magEarth(3)/magNormDip(idx,1));
end

eul6D = quat2Eul(output.quat6D);
eul9D = quat2Eul(output.quat9D);

plotMagNormDip = figure();
hold on; grid on;
xlabel('time_huro (s)'); ylabel('Magnetometer Norm Dip (deg)');
magNormDipPlot = plot(time_huro, magNormDip(:,2).*rad2deg, 'LineWidth', 1);



plotMagNorm = figure();
hold on; grid on;
xlabel('time_huro (s)'); ylabel('Norm Dip');
magNormDipPlot = plot(time_huro, magNormDip(:,1), 'LineWidth', 1);
%%
plotEul9D = figure();
hold on; grid on;
xlabel('time_huro (s)'); ylabel('Euler Angles (degrees)');
huro_rol = plot(time_huro, rawDataHuro.roll, 'r');
huro_pit = plot(time_huro, rawDataHuro.pitch, 'g');
huro_yaw = plot(time_huro, rawDataHuro.yaw, 'b');
vqf_rol = plot(time_huro, eul9D(:,1), 'r--');
vqf_pit = plot(time_huro, eul9D(:,2), 'g--');
vqf_yaw = plot(time_huro, eul9D(:,3), 'b--');
% xsen_rol = plot(time_xsen, rawDataXsen.roll, 'r--');
% xsen_pit = plot(time_xsen, rawDataXsen.pitch, 'g--');
% xsen_yaw = plot(time_xsen, rawDataXsen.yaw, 'b--');

legend([huro_rol, huro_pit, huro_yaw], {'Roll', 'Pitch', 'Yaw'});
title('Euler Angles from 9D Quaternion');

%%
function eul = quat2Eul(quat)
% Convert quaternion to Euler angles
rad2deg = 180 / pi;
q0 = quat(:, 1);
q1 = quat(:, 2);
q2 = quat(:, 3);
q3 = quat(:, 4);

roll = atan2(2*(q0.*q1 + q2.*q3), 1 - 2*(q1.^2 + q2.^2)) * rad2deg;
pitch = asin(2*(q0.*q2 - q3.*q1)) * rad2deg;
yaw = atan2(2*(q0.*q3 + q1.*q2), 1 - 2*(q2.^2 + q3.^2)) * rad2deg;

eul = [roll, pitch, yaw];
end

function out = quatRotate(q, v)
x = (1 - 2*q(3)^2 - 2*q(4)^2)*v(1) + 2*(q(3)*q(2) - q(1)*q(4))*v(2) + 2*(q(1)*q(3) + q(4)*q(2))*v(3);
y = 2*(q(1)*q(4) + q(3)*q(2))*v(1) + (1 - 2*q(2)^2 - 2*q(4)^2)*v(2) + 2*(q(3)*q(4) - q(2)*q(1))*v(3);
z = 2*(q(4)*q(2) - q(1)*q(3))*v(1) + 2*(q(1)*q(2) + q(4)*q(3))*v(2) + (1 - 2*q(2)^2 - 2*q(3)^2)*v(3);
out = [x; y; z];
end

