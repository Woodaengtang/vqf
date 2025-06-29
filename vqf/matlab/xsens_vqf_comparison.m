close all; clear all; clc;
slowRawData = readtable("log\comparison\XSens_06_29_1703.csv");
fastRawData = readtable("log\comparison\XSens_06_29_1704.csv");

ts = 1/100;

slow_vqf = VQF(ts, ts, ts);
fast_vqf = VQF(ts, ts, ts);

slow_raw_gry = [slowRawData.gx, slowRawData.gy, slowRawData.gz];
slow_raw_acc = [slowRawData.ax, slowRawData.ay, slowRawData.az];
slow_raw_mag = [slowRawData.mx, slowRawData.my, slowRawData.mz];

slow_output = slow_vqf.updateBatch(slow_raw_gry, slow_raw_acc, slow_raw_mag);
time_xsens = ts : ts : ts*length(slowRawData.roll);
slow_eul = quat2Eul(slow_output.quat6D);

plotSlowRol = figure();
hold on; grid on;
% vqf_rol = plot(time_xsens, slow_eul(:,1), "LineWidth",1);
% xsn_rol = plot(time_xsens, slowRawData.roll, "LineWidth",1);
% vqf_pit = plot(time_xsens, slow_eul(:,2), "LineWidth",1);
% xsn_pit = plot(time_xsens, slowRawData.pitch, "LineWidth",1);
vqf_yaw = plot(time_xsens, slow_eul(:,3), "LineWidth",1);
xsn_yaw = plot(time_xsens, slowRawData.yaw, "LineWidth",1);
% legend([vqf_rol, xsn_rol], {"VQF", "Xsens"});

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