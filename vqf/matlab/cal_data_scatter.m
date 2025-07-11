close all; clear all; clc;
debugData = readtable("log\Huro_07_11_1114.csv");
input_format = 'MM.dd HH:mm:ss:SSS';
debugData.Time = datetime(debugData.Time, 'InputFormat', input_format);

Scat = figure();
subplot(1, 2, 1);
hold on; grid on;
scatter3(debugData.debug0, debugData.debug1, debugData.debug2, 'Marker', '.');
axis equal;
subplot(1, 2, 2);
hold on; grid on;
scatter3(debugData.mx, debugData.my, debugData.mz, 'Marker', '.');
axis equal;
