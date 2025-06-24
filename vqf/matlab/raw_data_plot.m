close all; clear all; clc;

rawData = readtable("log\Huro_06_23_1636.csv");
diffData = readtable("log\Diff_06_23_1636.csv");
%%
time_step = 1/570;
data_len = length(rawData.Counter);
time = linspace(time_step, data_len*time_step, data_len);

rawAccPlot = figure();
hold on; grid on;
rawAx = plot(time, rawData.ax, 'LineWidth', 1);
rawAy = plot(time, rawData.ay, 'LineWidth', 1);
rawAz = plot(time, rawData.az, 'LineWidth', 1);
legend([rawAx, rawAy, rawAz], {"raw acc x", "raw acc y", "raw acc z"});
xlabel("time (s)"); ylabel("Acc (m/s^2)"); title("Raw data - Accelerometer");
print("raw_acc_plot_case_2.png", "-dpng", "-r500");

rad2deg = 180/pi;
rawGryPlot = figure();
hold on; grid on;
rawGx = plot(time, rawData.gx.*rad2deg, 'LineWidth', 1);
rawGy = plot(time, rawData.gy.*rad2deg, 'LineWidth', 1);
rawGz = plot(time, rawData.gz.*rad2deg, 'LineWidth', 1);
legend([rawGx, rawGy, rawGz], {"raw gry x", "raw gry y", "raw gry z"});
xlabel("time (s)"); ylabel("Gry (deg/s)"); title("Raw data - Gyroscope");
print("raw_gry_plot_case_2.png", "-dpng", "-r500");

rawMagPlot = figure();
hold on; grid on;
rawMx = plot(time, rawData.mx, 'LineWidth', 1);
rawMy = plot(time, rawData.my, 'LineWidth', 1);
rawMz = plot(time, rawData.mz, 'LineWidth', 1);
rawMnorm = plot(time, vecnorm([rawData.mx, rawData.my, rawData.mz], 2, 2));
legend([rawMx, rawMy, rawMz, rawMnorm], {"raw mag x", "raw mag y", "raw mag z", "raw mag norm"}, "Location", "best");
xlabel("time (s)"); ylabel("Mag (uT)"); title("Raw data - Magnetometer");
print("raw_mag_plot_case_2.png", "-dpng", "-r500");

rawMagXY = figure();
hold on; grid on; axis equal;
magXY = plot(rawData.mx, rawData.my, 'LineWidth', 1);
estK = scatter(mean(rawData.mx), mean(rawData.my), 'r', 'filled');
xline(0, 'k--'); yline(0, 'k--');
legend([magXY, estK], {"raw mag xy", "estimated position of rotation axis"});
xlabel("Mag x (uT)"); ylabel("Mag y (uT)"); title("Raw data - Magnetometer XY plane");
print("raw_mag_xy_plot_case_2.png", "-dpng", "-r500");

vqfRPY = figure();
hold on; grid on;
rawR = plot(time, rawData.roll, 'LineWidth', 1);
rawP = plot(time, rawData.pitch, 'LineWidth', 1);
rawY = plot(time, rawData.yaw, 'LineWidth', 1);
legend([rawR, rawP, rawY], {"raw roll", "raw pitch", "raw yaw"});
xlabel("time (s)"); ylabel("Angle (degrees)"); title("Raw data - Roll, Pitch, Yaw");
print("raw_rpy_plot_case_2.png", "-dpng", "-r500");

%%
rawX = rawData.mx;
rawY = rawData.my;

y = -rawX.^2;
c1 = rawX.*rawY;
c2 = rawY.^2;
c3 = rawX;
c4 = rawY;
c5 = ones([data_len, 1]);

A = [c1, c2, c3, c4, c5];
coeff = ((A'*A)\A') * y;

a_ = 1;
b_ = coeff(1);
c_ = coeff(2);
d_ = coeff(3);
e_ = coeff(4);
f_ = coeff(5);

delta = atan2(b_, 1 - c_) / 2;
c_x = (2*c_*d_ - b_*e_) / (b_^2 - 4 * c_);
c_y = (2*a_*e_ - b_*d_) / (b_^2 - 4 * c_);
cos_d = cos(delta);
sin_d = sin(delta);
w = sqrt((c_x^2 + b_*c_x*c_y  + c_*c_y^2 - f_) / (cos_d^2 + b_*cos_d*sin_d + c_*sin_d^2));
h = sqrt((c_x^2 + b_*c_x*c_y  + c_*c_y^2 - f_) / (sin_d^2 + b_*cos_d*sin_d + c_*cos_d^2));
sigma = h/w;
tmpxcal = (cos_d * (rawX - c_x) + sin_d * (rawY - c_y)) * sigma;
tmpycal = -sin_d * (rawX - c_x) + cos_d * (rawY - c_y);
xcal = cos_d * tmpxcal - sin_d * tmpycal;
ycal = sin_d * tmpxcal + cos_d * tmpycal;

calMagXY = figure();
hold on; 
rawXY = scatter(rawData.mx, rawData.my, 'Marker', '.');
calXY = scatter(xcal, ycal, 'Marker', '.');
xline(0, 'k--'); yline(0, 'k--');
grid on; 
legend([rawXY, calXY], {"raw mag xy", "cal mag xy"});
axis equal;
xlabel("Mag x (uT)"); ylabel("Mag y (uT)"); title("comparison of raw and calibrated magnetometer data in XY plane");
print("raw_cal_mag_xy.png", "-dpng", "-r500");

ax2_o = -rawData.mx.^2;
ay2_i = rawData.my.^2;
az2_i = rawData.mz.^2;
axy_i = rawData.mx.*rawData.my;
ayz_i = rawData.my.*rawData.mz;
azx_i = rawData.mz.*rawData.mx;
ax_i  = rawData.mx;
ay_i  = rawData.my;
az_i  = rawData.mz;
a0_i  = ones([length(rawData.mx), 1]);

A3d = [ay2_i, az2_i, axy_i, ayz_i, azx_i, ax_i, ay_i, az_i, a0_i];

coeff3d = (A3d'*A3d)\A3d'*ax2_o;

center = [         -2,   -coeff3d(3),   -coeff3d(5);...
 -coeff3d(3), -2*coeff3d(1),   -coeff3d(4);...
 -coeff3d(5),   -coeff3d(4), -2*coeff3d(2)]\[coeff3d(6); coeff3d(7); coeff3d(8)];

% xyplane.A = []

% A_ = [ax_i-center(1), ay_i-center(2), az_i-center(3)];
% [V_, D_] = eig(A_'*A_);
% 
% calA = A_*V_;

r = 35;
[sp_x, sp_y, sp_z] = sphere(100);
sp_x = sp_x*r;
sp_y = sp_y*r;
sp_z = sp_z*r;

% rawMagXYZ = figure();
% hold on; grid on;
% % rawXYZ = scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
% calXYZ = scatter3(ax_i-center(1), ay_i-center(2), az_i-center(3), 'Marker', '.');
% surf(sp_x, sp_y, sp_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.1);
% axis equal; 
% xlabel('x'); ylabel('y'); zlabel('z');
% view([-51, 32]);


% rawGyro = figure();
% hold on; grid on;
% rawGryx = plot(time, rawData.gx, 'r');
% rawGryy = plot(time, rawData.gy, 'g');
% rawGryz = plot(time, rawData.gz, 'b');


buf_hz = 1/0.05;
buf_time = 1/buf_hz : 1/buf_hz : 1/buf_hz*length(diffData.Time);
diffPlot = figure();
hold on; grid on;
diffRol = plot(buf_time, diffData.Diff_Roll, 'r', 'LineWidth', 1);
diffPit = plot(buf_time, diffData.Diff_Pitch, 'g', 'LineWidth', 1);
diffYaw = plot(buf_time, diffData.Diff_Yaw, 'b', 'LineWidth', 1);
ylim([-6, 6]);

% huroRol = plot(buf_time, diffData.Huro_Roll, 'r');
% huroPit = plot(buf_time, diffData.Huro_Pitch, 'g');
% huroYaw = plot(buf_time, diffData.Huro_Yaw, 'b');
% xsenRol = plot(buf_time, diffData.XSens_Roll, 'r--');
% xsenPit = plot(buf_time, diffData.XSens_Pitch, 'g--');
% xsenYaw = plot(buf_time, diffData.XSens_Yaw, 'b--');

% 
% compPlot = figure();
% hold on; grid on;
% huroRol = plot(buf_time, diffData.Huro_Roll, 'r');
% huroPit = plot(buf_time, diffData.Huro_Pitch, 'g');
% huroYaw = plot(buf_time, diffData.Huro_Yaw, 'b');
% xsenRol = plot(buf_time, diffData.XSens_Roll, 'r--');
% xsenPit = plot(buf_time, diffData.XSens_Pitch, 'g--');
% xsenYaw = plot(buf_time, diffData.XSens_Yaw, 'b--');