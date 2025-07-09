close all; clear all; clc;
rawData = readtable("log\Huro_06_20_1548.csv");

time_step = 1/570;
data_len = length(rawData.Counter);
time = linspace(time_step, data_len*time_step, data_len);

y = -rawData.mx.^2;
c1 = rawData.mx.*rawData.my;
c2 = rawData.my.^2;
c3 = rawData.mx;
c4 = rawData.my;
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
tmpxcal = (cos_d * (rawData.mx - c_x) + sin_d * (rawData.my - c_y)) * sigma;
tmpycal = -sin_d * (rawData.mx - c_x) + cos_d * (rawData.my - c_y);
xcal = cos_d * tmpxcal - sin_d * tmpycal;
ycal = sin_d * tmpxcal + cos_d * tmpycal;

calMagXY = figure();
hold on; 
rawDataXY = scatter(rawData.mx, rawData.my, 'Marker', '.');
calXY = scatter(xcal, ycal, 'Marker', '.');
xline(0, 'k--'); yline(0, 'k--');
grid on; 
legend([rawData.mxY, calXY], {"raw mag xy", "cal mag xy"});
axis equal;
xlabel("Mag x (uT)"); ylabel("Mag y (uT)"); title("comparison of raw and calibrated magnetometer data in XY plane");
print("raw_cal_mag_xy.png", "-dpng", "-r500");
