close all; clear all; clc;

% rawData = readtable("Huro_06_17_1059.csv");
rawData = readtable("Huro_06_17_1316.csv");

time_step = 0.01;
data_len = length(rawData.Counter);
time = linspace(time_step, data_len*time_step, data_len);

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

xcal = cos_d * (rawX - c_x) + sin_d * (rawY - c_y);
ycal = -sin_d * (rawX - c_x) + cos_d * (rawY - c_y);
xcal = xcal / w * h;

