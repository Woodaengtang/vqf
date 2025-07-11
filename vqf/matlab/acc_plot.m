close all; clear all; clc;
readData = readtable("log\Huro_07_10_1632.csv");
input_format = 'MM.dd HH:mm:ss:SSS';
readData.Time = datetime(readData.Time, 'InputFormat', input_format);

figParam.line_width = 1.5;

accPlot = figure();
subplot(3, 1, 1);
hold on; grid on;
plot(readData.Time, readData.ax, 'LineWidth', figParam.line_width);
subplot(3, 1, 2);
hold on; grid on;
plot(readData.Time, readData.ay, 'LineWidth', figParam.line_width);
subplot(3, 1, 3);
hold on; grid on;
plot(readData.Time, readData.az, 'LineWidth', figParam.line_width);

gyrPlot = figure();
subplot(3, 1, 1);
hold on; grid on;
plot(readData.Time, readData.gx, 'LineWidth', figParam.line_width);
subplot(3, 1, 2);
hold on; grid on;
plot(readData.Time, readData.gy, 'LineWidth', figParam.line_width);
subplot(3, 1, 3);
hold on; grid on;
plot(readData.Time, readData.gz, 'LineWidth', figParam.line_width);