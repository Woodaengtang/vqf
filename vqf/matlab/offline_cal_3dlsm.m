close all; clear all; clc;
% rawData = readtable("log\cal_data\Huro_07_04_1626.csv");
% rawData = readtable("log\Huro_07_08_1455.csv");
rawData = readtable("log\Huro_07_10_1654.csv");
input_format = 'MM.dd HH:mm:ss:SSS';
rawData.Time = datetime(rawData.Time, 'InputFormat', input_format);

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
y = ones(height(rawData), 1);

A = [mc_1 mc_2 mc_3 mc_4 mc_5 mc_6 mc_7 mc_8 mc_9];
x = (A' * A) \ (A' * y);
M = [x(1) x(4) x(5);
    x(4) x(2) x(6);
    x(5) x(6) x(3)];
b = [x(7), x(8), x(9)];

[V, D] = eig(M);

new_b = b * V;

% Scale the matrix to make it a unit sphere
scale_q = sqrt(new_b(1)^2/(4*D(1, 1)) + new_b(2)^2/(4*D(2, 2)) + new_b(3)^2/(4*D(3, 3)) + 1);
scale_matrix = sqrt(D)./scale_q;

hard_iron = [new_b(1)/D(1,1), new_b(2)/D(2, 2), new_b(3)/D(3, 3)]./2;
soft_iron = V*scale_matrix;
cal_mag_xyz = (soft_iron * (V' * [rawData.mx, rawData.my, rawData.mz]' + hard_iron'))';

str_V = sprintf('%.4ff, ', V'); 
str_V = str_V(1:end-2);
fprintf('vqf_real_t V[%d] = {%s};\n', numel(V), str_V);
str_hard_iron = sprintf('%.4ff, ', hard_iron);
str_hard_iron = str_hard_iron(1:end-2);
fprintf('vqf_real_t hard_iron[%d] = {%s};\n', numel(hard_iron), str_hard_iron);
str_soft_iron = sprintf('%.4ff, ', soft_iron');
str_soft_iron = str_soft_iron(1:end-2);
fprintf('vqf_real_t soft_iron[%d] = {%s};\n', numel(soft_iron), str_soft_iron);

view_angle = [-30, 20];
calPlot = figure();
% calPlot.Position = 1.0e+03.*[0.0010    0.0490    1.9200    0.9568];
subplot(2, 4, 1);
grid on; hold on; 
scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
xlabel('mx'); ylabel('my'); zlabel('mz');
axis equal;
view(view_angle);
subplot(2, 4, 2);
grid on; hold on; 
scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
xlabel('mx'); ylabel('my'); zlabel('mz');
axis equal;
view([0, 90]);
subplot(2, 4, 5);
grid on; hold on; 
scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
xlabel('mx'); ylabel('my'); zlabel('mz');
axis equal;
view([90, 0]);
subplot(2, 4, 6);
grid on; hold on; 
scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
xlabel('mx'); ylabel('my'); zlabel('mz');
axis equal;
view([0, 0]);
subplot(2, 4, [3, 4, 7, 8]);
grid on; hold on;
calMag = scatter3(cal_mag_xyz(:, 1), cal_mag_xyz(:, 2), cal_mag_xyz(:, 3), 'Marker', '.', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
xlabel('mx'); ylabel('my'); zlabel('mz');
subtitle('Calibrated Mag');
axis equal;
view(view_angle);


datalength = length(rawData.Time)-1;
if rem(datalength, 2) == 0
    h_len = datalength/2;
else
    h_len = (datalength-1)/2;
end

idx = min(find(rawData.roll, 1, 'first'));
az = 32;
el = 30;
plotDiff = figure();
plotDiff.Position = [488, 132.2, 743.4, 529.8];
hold on; grid on;
scatter3(rawData.debug0, rawData.debug1, rawData.debug2, 'Marker', '.');
% raw = scatter3(rawData.mx(1:idx), rawData.my(1:idx), rawData.mz(1:idx), 'Marker', '.');
% cal = scatter3(rawData.mx(idx+1:end), rawData.my(idx+1:end), rawData.mz(idx+1:end), 'Marker', '.');
% xlabel('mx'); ylabel('my'); zlabel('mz'); legend([raw, cal], {'raw data', 'calibrated data'});
title('Calibration process validation');
view([az, el]);

