close all; clear all; clc;
rawData = readtable("log\Huro_06_30_1418.csv");

mean_gx = mean(rawData.gx);
mean_gy = mean(rawData.gy);
mean_gz = mean(rawData.gz);


ts = 1/570;
time = linspace(ts, ts*length(rawData.gx), length(rawData.gx));

rawGryPlot = figure();
rawGryPlot.Position = [488.0000, 49.8000, 561.0000, 948.8000];
subplot(3, 1, 1);
hold on; grid on;
rawGx  = plot(time, rawData.gx);
% meanGx = plot(time, mean_gx*ones([length(time), 1]), 'r--', 'LineWidth', 2);
estBiasGx = plot(time, rawData.debug0, 'g--', 'LineWidth', 1);
rawData.debug0
xlabel('time (s)'); ylabel('gx');
subplot(3, 1, 2);
hold on; grid on;
rawGy  = plot(time, rawData.gy);
% meanGy = plot(time, mean_gy*ones([length(time), 1]), 'r--', 'LineWidth', 2);
estBiasGy = plot(time, rawData.debug1, 'g--', 'LineWidth', 1);
xlabel('time (s)'); ylabel('gy');
subplot(3, 1, 3);
hold on; grid on;
rawGz  = plot(time, rawData.gz);
% meanGz = plot(time, mean_gz*ones([length(time), 1]), 'r--', 'LineWidth', 2);
estBiasGz = plot(time, rawData.debug2, 'g--', 'LineWidth', 1);
xlabel('time (s)'); ylabel('gz');


fprintf('vqf_real_t init_cal[3] = {%.6ff, %.6ff, %.6ff};\n', mean_gx, mean_gy, mean_gz);