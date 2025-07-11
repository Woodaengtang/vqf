close all; clear all; clc;
%% Base data
% diffData = readtable("log\Diff_07_10_1005.csv");
% diffData = readtable("log\Diff_07_10_1042.csv");
% diffData = readtable("log\Diff_07_10_1104.csv");
% diffData = readtable("log\Diff_07_10_1105.csv");
% diffData = readtable("log\Diff_07_10_1120.csv");
%% Cali, delta data
% diffData = readtable("log\Diff_07_10_1150.csv");
% diffData = readtable("log\Diff_07_10_1315.csv");
% diffData = readtable("log\Diff_07_10_1326.csv");
% diffData = readtable("log\Diff_07_10_1335.csv");
diffData = readtable("log\Diff_07_10_1632.csv");

input_format = 'MM.dd HH:mm:ss:SSS';
diffData.Time = datetime(diffData.Time, 'InputFormat', input_format);

for i = 1:length(diffData.Time)
    if diffData.Diff_Yaw(i) > 180
        diffData.Diff_Yaw(i) = diffData.Diff_Yaw(i) - 360;
    elseif diffData.Diff_Yaw(i) < -180
        diffData.Diff_Yaw(i) = diffData.Diff_Yaw(i) + 360;
    end
end

figParam.angle_lim = 10;
figParam.line_width = 1.5;
x_coord = [diffData.Time(1), diffData.Time(end)];
x_patch = [x_coord(1), x_coord(2), x_coord(2), x_coord(1)];
y_coord = [-2, 2];
y_patch = [y_coord(1), y_coord(1), y_coord(2), y_coord(2)];

diffAnglePlot = figure();
diffAnglePlot.Position = [108.2, 913.8, 862.4, 948.8];
subplot(3, 1, 1);
title('Angle difference btw Xsens and Huro');
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(diffData.Time, diffData.Diff_Roll, 'LineWidth', figParam.line_width);
xlim([diffData.Time(1), diffData.Time(end)]); ylim([-figParam.angle_lim, figParam.angle_lim]);
ylabel('Roll (deg)');
subplot(3, 1, 2);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(diffData.Time, diffData.Diff_Pitch, 'LineWidth', figParam.line_width);
xlim([diffData.Time(1), diffData.Time(end)]); ylim([-figParam.angle_lim, figParam.angle_lim]);
ylabel('Pitch (deg)');
subplot(3, 1, 3);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(diffData.Time, diffData.Diff_Yaw, 'LineWidth', figParam.line_width);
xlim([diffData.Time(1), diffData.Time(end)]); 
% ylim([-figParam.angle_lim, figParam.angle_lim]);
ylabel('Yaw (deg)');

diffXsensHuro = figure();
diffXsensHuro.Position = [5.0000, 913.8000, 856.0000, 948.8000];
subplot(2, 1, 1);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(diffData.Time, diffData.Diff_Yaw, 'LineWidth', figParam.line_width);
xlim([diffData.Time(1), diffData.Time(end)]); ylabel('Yaw (deg)');
subplot(2, 1, 2);
hold on; grid on;
yaw_huro = plot(diffData.Time, diffData.Huro_Yaw, 'LineWidth', figParam.line_width);
yaw_xsens = plot(diffData.Time, diffData.XSens_Yaw, 'LineWidth', figParam.line_width);
legend([yaw_huro, yaw_xsens], {'\psi Huro', '\psi Xsens'}); ylabel('Yaw (deg)');

heading_rmse = rmse(zeros([length(diffData.Diff_Yaw), 1]), diffData.Diff_Yaw);
disp(heading_rmse);

