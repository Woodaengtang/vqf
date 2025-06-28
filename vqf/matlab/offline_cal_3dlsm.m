close all; clear all; clc;
rawData = readtable("log\Huro_06_23_1544.csv");

% mc : matrix column
mc_1 = rawData.mx.^2;
mc_2 = rawData.my.^2;
mc_3 = rawData.mz.^2;
mc_4 = 2*rawData.mx.*rawData.my;
mc_5 = 2*rawData.mx.*rawData.mz;
mc_6 = 2*rawData.my.*rawData.mz;
mc_7 = rawData.mx;
mc_8 = rawData.my;
mc_9 = rawData.mz;

% y : output 
y = -1*ones(height(rawData), 1);

A = [mc_1 mc_2 mc_3 mc_4 mc_5 mc_6 mc_7 mc_8 mc_9];
x = (A' * A) \ (A' * y);
M = [x(1) x(4) x(5);
    x(4) x(2) x(6);
    x(5) x(6) x(3)];
b = [x(7), x(8), x(9)];

[V, D] = eig(M);

new_b = b * V;

% Scale the matrix to make it a unit sphere
scale_q = sqrt(new_b(1)^2/(4*D(1, 1)) + new_b(2)^2/(4*D(2, 2)) + new_b(3)^2/(4*D(3, 3)) - 1);
scale_matrix = sqrt(D)./scale_q;

hard_iron = [new_b(1)/D(1,1), new_b(2)/D(2, 2), new_b(3)/D(3, 3)]./2;
soft_iron = V*scale_matrix;
cal_mag_xyz = (soft_iron * (V' * [rawData.mx, rawData.my, rawData.mz]' + hard_iron'))'./scale_matrix(1,1);


rawPlot = figure();
grid on; hold on; 
rawMag = scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
calMag = scatter3(cal_mag_xyz(:, 1), cal_mag_xyz(:, 2), cal_mag_xyz(:, 3), 'Marker', '.');
axis equal;


