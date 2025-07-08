close all; clear all; clc;
rawData = readtable("log\Huro_07_08_1555.csv");
input_format = 'MM.dd HH:mm:ss:SSS';
rawData.Time = datetime(rawData.Time, 'InputFormat', input_format);


angle = atan2(rawData.debug0, rawData.debug1);
new_angle = NaN([length(angle), 1]);
new_angle(1) = angle(1);
alpha = 0.9512;
for idx = 1:length(angle)-1
    new_angle(idx+1) = (1-alpha)*angle(idx+1) + alpha*new_angle(idx);
end
rad2deg_c = 180/pi;
magDebug = figure();
magDebug.Position = 1.0e+03 .*[0.4650, 1.1114, 0.8416, 0.4200];
hold on; grid on;
idx = 10997;
raw_heading = plot(rawData.Time(idx:end), angle(idx:end).*rad2deg_c);
lpf = plot(rawData.Time(idx:end), new_angle(idx:end).*rad2deg_c, 'LineWidth', 1);
delta = plot(rawData.Time(idx:end), rawData.debug3(idx:end).*rad2deg_c, 'LineWidth', 1);
ylabel('heading angle (deg)'); 
legend([raw_heading, lpf, delta], {"raw heading \psi", "LPF(\psi) 5Hz", "\delta"});

%%
