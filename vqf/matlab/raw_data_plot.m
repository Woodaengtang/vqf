close all; clear all; clc;

rawData = readtable("log\Huro_07_10_1005.csv");
xsenData = readtable("log\XSens_06_29_1126.csv");
%%
time_step = 1/570;
data_len = length(rawData.Counter);
time = linspace(time_step, data_len*time_step, data_len);
x_time_step = 1/100;
x_data_len = length(xsenData.mx);
x_time = linspace(x_time_step, x_data_len*x_time_step, x_data_len);
% rawAccPlot = figure();
% hold on; grid on;
% rawAx = plot(time, rawData.ax, 'LineWidth', 1);
% rawAy = plot(time, rawData.ay, 'LineWidth', 1);
% rawAz = plot(time, rawData.az, 'LineWidth', 1);
% legend([rawAx, rawAy, rawAz], {"raw acc x", "raw acc y", "raw acc z"});
% xlabel("time (s)"); ylabel("Acc (m/s^2)"); title("Raw data - Accelerometer");
% print("raw_acc_plot_case_2.png", "-dpng", "-r500");
% 
% rad2deg = 180/pi;
% rawGryPlot = figure();
% hold on; grid on;
% rawGx = plot(time, rawData.gx.*rad2deg, 'LineWidth', 1);
% rawGy = plot(time, rawData.gy.*rad2deg, 'LineWidth', 1);
% rawGz = plot(time, rawData.gz.*rad2deg, 'LineWidth', 1);
% legend([rawGx, rawGy, rawGz], {"raw gry x", "raw gry y", "raw gry z"});
% xlabel("time (s)"); ylabel("Gry (deg/s)"); title("Raw data - Gyroscope");
% print("raw_gry_plot_case_2.png", "-dpng", "-r500");

rawMagPlot = figure();
hold on; grid on;
rawMx = plot(time, rawData.mx, 'LineWidth', 1);
% xsnMx = plot(x_time, xsenData.mx, 'LineWidth', 1);
rawMy = plot(time, rawData.my, 'LineWidth', 1);
rawMz = plot(time, rawData.mz, 'LineWidth', 1);
rawMnorm = plot(time, vecnorm([rawData.mx, rawData.my, rawData.mz], 2, 2));
legend([rawMx, rawMy, rawMz, rawMnorm], {"raw mag x", "raw mag y", "raw mag z", "raw mag norm"}, "Location", "best");
xlabel("time (s)"); ylabel("Mag (uT)"); title("Raw data - Magnetometer");
print("raw_mag_plot_case_2.png", "-dpng", "-r500");

% rawMagXY = figure();
% hold on; grid on; axis equal;
% magXY = plot(rawData.mx, rawData.my, 'LineWidth', 1);
% estK = scatter(mean(rawData.mx), mean(rawData.my), 'r', 'filled');
% xline(0, 'k--'); yline(0, 'k--');
% legend([magXY, estK], {"raw mag xy", "estimated position of rotation axis"});
% xlabel("Mag x (uT)"); ylabel("Mag y (uT)"); title("Raw data - Magnetometer XY plane");
% print("raw_mag_xy_plot_case_2.png", "-dpng", "-r500");
% 
% vqfRPY = figure();
% hold on; grid on;
% rawR = plot(time, rawData.roll, 'LineWidth', 1);
% rawP = plot(time, rawData.pitch, 'LineWidth', 1);
% rawY = plot(time, rawData.yaw, 'LineWidth', 1);
% legend([rawR, rawP, rawY], {"raw roll", "raw pitch", "raw yaw"});
% xlabel("time (s)"); ylabel("Angle (degrees)"); title("Raw data - Roll, Pitch, Yaw");
% print("raw_rpy_plot_case_2.png", "-dpng", "-r500");
% 
% 
% buf_hz = 1/0.05;
% buf_time = 1/buf_hz : 1/buf_hz : 1/buf_hz*length(diffData.Time);
% diffPlot = figure();
% hold on; grid on;
% diffRol = plot(buf_time, diffData.Diff_Roll, 'r', 'LineWidth', 1);
% diffPit = plot(buf_time, diffData.Diff_Pitch, 'g', 'LineWidth', 1);
% diffYaw = plot(buf_time, diffData.Diff_Yaw, 'b', 'LineWidth', 1);
% ylim([-6, 6]);
% 
% huroRol = plot(buf_time, diffData.Huro_Roll, 'r');
% huroPit = plot(buf_time, diffData.Huro_Pitch, 'g');
% huroYaw = plot(buf_time, diffData.Huro_Yaw, 'b');
% xsenRol = plot(buf_time, diffData.XSens_Roll, 'r--');
% xsenPit = plot(buf_time, diffData.XSens_Pitch, 'g--');
% xsenYaw = plot(buf_time, diffData.XSens_Yaw, 'b--');
% 
% 
% compPlot = figure();
% hold on; grid on;
% huroRol = plot(buf_time, diffData.Huro_Roll, 'r');
% huroPit = plot(buf_time, diffData.Huro_Pitch, 'g');
% huroYaw = plot(buf_time, diffData.Huro_Yaw, 'b');
% xsenRol = plot(buf_time, diffData.XSens_Roll, 'r--');
% xsenPit = plot(buf_time, diffData.XSens_Pitch, 'g--');
% xsenYaw = plot(buf_time, diffData.XSens_Yaw, 'b--');