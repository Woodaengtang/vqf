close all; clear all; clc;
% rawData = readtable("log\Huro_07_11_1348.csv");
% rawData = readtable("log\Huro_07_11_1356.csv");
% rawData = readtable("log\Huro_07_11_1404.csv");
rawData = readtable("log\Huro_07_11_1414.csv");
% rawData = readtable("log\Huro_07_11_1421.csv");
% rawData = readtable("log\Huro_07_11_1430.csv");
% rawData = readtable("log\Huro_07_10_1315.csv");

input_format = 'MM.dd HH:mm:ss:SSS';
rawData.Time = datetime(rawData.Time, 'InputFormat', input_format);

cutoff_freq = 5;
sampling_time = 1/570;
alpha = lpf_mag(cutoff_freq, sampling_time);

new_mag = NaN([length(rawData.mx), 3]);
new_mag(1,:) = [rawData.mx(1), rawData.mx(2), rawData.mx(3)];
for idx = 2:length(new_mag)
    new_mag(idx, 1) = (1-alpha)*rawData.mx(idx) + alpha*new_mag(idx-1, 1);
    new_mag(idx, 2) = (1-alpha)*rawData.my(idx) + alpha*new_mag(idx-1, 2);
    new_mag(idx, 3) = (1-alpha)*rawData.mz(idx) + alpha*new_mag(idx-1, 3);
end

angle = atan2(rawData.debug0, rawData.debug1);
new_angle = NaN([length(angle), 1]);
new_angle(1) = angle(1);

for idx = 1:length(angle)-1
    new_angle(idx+1) = (1-alpha)*angle(idx+1) + alpha*new_angle(idx);
end
rad2deg_c = 180/pi;
% 
magDebug = figure();
magDebug.Position = 1.0e+03 .*[0.4650, 1.1114, 0.8416, 0.4200];
hold on; grid on;
idx = 10997;
% raw_heading = plot(rawData.Time(idx:end), angle(idx:end).*rad2deg_c);
% delta = plot(rawData.Time(idx:end), rawData.debug3(idx:end).*rad2deg_c, 'LineWidth', 1.5);
raw_heading = plot(rawData.Time, angle.*rad2deg_c);
delta = plot(rawData.Time, rawData.debug3.*rad2deg_c, 'LineWidth', 1.5);
ylabel('heading angle (deg)'); 
legend([raw_heading, delta], {"raw heading \psi", "\delta"});

lpf_mag_fig = figure();
subplot(3, 1, 1);
hold on; grid on;
% plot(rawData.Time(idx:end), rawData.mx(idx:end));
% plot(rawData.Time(idx:end), new_mag(idx:end, 1), 'LineWidth', 1);
plot(rawData.Time, rawData.mx);
plot(rawData.Time, new_mag(:, 1), 'LineWidth', 1);
subplot(3, 1, 2);
hold on; grid on;
% plot(rawData.Time(idx:end), rawData.my(idx:end));
% plot(rawData.Time(idx:end), new_mag(idx:end, 2), 'LineWidth', 1);
plot(rawData.Time, rawData.my);
plot(rawData.Time, new_mag(:, 2), 'LineWidth', 1);
subplot(3, 1, 3);
hold on; grid on;
% plot(rawData.Time(idx:end), rawData.mz(idx:end));
% plot(rawData.Time(idx:end), new_mag(idx:end, 3), 'LineWidth', 1);
plot(rawData.Time, rawData.mz);
plot(rawData.Time, new_mag(:, 3), 'LineWidth', 1);


%%
function alpha = lpf_mag(cutoff_freq, sampling_time)
    alpha = exp(-cutoff_freq*sampling_time*2*pi);
end