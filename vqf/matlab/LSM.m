function param = mag_cali_lsm(mag_x, mag_y)
%--------------------------------------------------------------------------
% LSM - Least Squares Method for magnetometer data fitting
%
% This function performs ellipse fitting using the Least Squares Method (LSM).
% It is designed to fit the general conic section equation:
%     ax^2 + bxy + cy^2 + dx + ey + f = 0
% to a set of 2D magnetometer data points.
%
% INPUT:
%   mag_x : Vector containing the x-axis data from the magnetometer.
%   mag_y : Vector containing the y-axis data from the magnetometer.
%
% The function estimates the parameters [b, c, d, e, f]/a.
%
%--------------------------------------------------------------------------
output = -mag_x.^2;
column1 = mag_x .* mag_y;
column2 = mag_y.^2;
column3 = mag_x;
column4 = mag_y;
column5 = ones(size(mag_x));
A = [column1, column2, column3, column4, column5];
coeff = (A'*A) \ A' * output;

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

param.c_x = c_x;
param.c_y = c_y;
param.cos_d = cos_d;
param.sin_d = sin_d;
param.sigma = h/w;

end