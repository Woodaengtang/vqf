close all; clear all; clc;
rand1 = readtable("log\cal_data\Huro_06_30_0939.csv");
rand2 = readtable("log\cal_data\Huro_06_30_1212.csv");
rand3 = readtable("log\cal_data\Huro_06_30_1545.csv");

view_point = [36, 36];
randPlot = figure();
hold on; grid on;
scatter3(rand1.mx, rand1.my, rand1.mz, '.');
scatter3(rand2.mx, rand2.my, rand2.mz, '.');
scatter3(rand3.mx, rand3.my, rand3.mz, '.');
axis equal;
view(view_point);
xlabel("mx"); ylabel("my"); zlabel("mz");