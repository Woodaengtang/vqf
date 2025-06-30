close all; clear all; clc;
rawData = readtable("log\mag_biases\Huro_06_29_1533.csv");
% rawData = readtable("log\mag_biases\Huro_06_29_1535.csv");
% rawData = readtable("log\mag_biases\Huro_06_29_1536.csv");

rest_idx = rawData.Rest(1);
rest_idx_prev = rest_idx;
count = 0;
ts = 1/570;
time = linspace(ts, ts*length(rawData.mx), length(rawData.mx));

moveidx = NaN([2, 1]);
for i = 1:length(rawData.Rest)
    rest_idx = rawData.Rest(i);
    if ~(rest_idx == rest_idx_prev)
        moveidx(count + 1) = i;
        if count == 0
            init_mx_mean = mean(rawData.mx(1:i-1));
            init_my_mean = mean(rawData.my(1:i-1));
            init_mz_mean = mean(rawData.mz(1:i-1));
            init_mx_std = std(rawData.mx(1:i-1));
            init_my_std = std(rawData.my(1:i-1));
            init_mz_std = std(rawData.mz(1:i-1));
        else
            mv_mx_mean = mean(rawData.mx(i:end));
            mv_my_mean = mean(rawData.my(i:end));
            mv_mz_mean = mean(rawData.mz(i:end));
            mv_mx_std = std(rawData.mx(i:end));
            mv_my_std = std(rawData.my(i:end));
            mv_mz_std = std(rawData.mz(i:end));
        end
        count = count + 1;
    end
    rest_idx_prev = rest_idx;
end

biasCheck = figure();
subplot([3, 1, 1]);
grid on; hold on;
plot(time, rawData.mx, 'LineWidth', 1);

subplot([3, 1, 2]);
grid on; hold on;
plot(time, rawData.my, 'LineWidth', 1);
subplot([3, 1, 3]);
plot(time, rawData.mz, 'LineWidth', 1);

init_mag = [init_mx_mean, init_my_mean, init_mz_mean, init_mx_std, init_my_std, init_mz_std];
mv_mag = [mv_mx_mean, mv_my_mean, mv_mz_mean, mv_mx_std, mv_my_std, mv_mz_std];
