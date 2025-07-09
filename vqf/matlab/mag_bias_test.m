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
        moveidx(count + 1) = i-1;
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
biasCheck.Position = [1000, 99, 839, 1139];
subplot(3, 1, 1);
grid on; hold on;
y_lims_init = [init_mx_mean-init_mx_std*3, init_mx_mean+init_mx_std*3];
y_coords_init = [y_lims_init(1), y_lims_init(1), y_lims_init(2), y_lims_init(2)];
x_patch1 = [time(1), time(moveidx(1)), time(moveidx(1)), time(1)];
patch(x_patch1, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
y_lims_mv = [mv_mx_mean-mv_mx_std*3, mv_mx_mean+mv_mx_std*3];
y_coords_init = [y_lims_mv(1), y_lims_mv(1), y_lims_mv(2), y_lims_mv(2)];
x_patch2 = [time(moveidx(2)), time(end), time(end), time(moveidx(2))];
patch(x_patch2, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time, rawData.mx);
plot(time(1:moveidx(1)), init_mx_mean*ones([moveidx(1), 1]), 'r--', 'LineWidth', 2);
plot(time(moveidx(2):end), mv_mx_mean*ones([length(time)-moveidx(2)+1, 1]), 'r--', 'LineWidth', 2);
y_lims = ylim;
y_coords= [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];
x_patch = [time(moveidx(1)), time(moveidx(2)), time(moveidx(2)), time(moveidx(1))];
patch(x_patch, y_coords, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('time (s)'); ylabel("mx");

subplot(3, 1, 2);
grid on; hold on;
y_lims_init = [init_my_mean-init_my_std*3, init_my_mean+init_my_std*3];
y_coords_init = [y_lims_init(1), y_lims_init(1), y_lims_init(2), y_lims_init(2)];
x_patch1 = [time(1), time(moveidx(1)), time(moveidx(1)), time(1)];
patch(x_patch1, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
y_lims_mv = [mv_my_mean-mv_my_std*3, mv_my_mean+mv_my_std*3];
y_coords_init = [y_lims_mv(1), y_lims_mv(1), y_lims_mv(2), y_lims_mv(2)];
x_patch2 = [time(moveidx(2)), time(end), time(end), time(moveidx(2))];
patch(x_patch2, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time, rawData.my);
plot(time(1:moveidx(1)), init_my_mean*ones([moveidx(1), 1]), 'r--', 'LineWidth', 1);
plot(time(moveidx(2):end), mv_my_mean*ones([length(time)-moveidx(2)+1, 1]), 'r--', 'LineWidth', 2);
y_lims = ylim;
y_coords= [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];
x_patch = [time(moveidx(1)), time(moveidx(2)), time(moveidx(2)), time(moveidx(1))];
patch(x_patch, y_coords, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('time (s)'); ylabel("my");

subplot(3, 1, 3);
grid on; hold on;
y_lims_init = [init_mz_mean-init_mz_std*3, init_mz_mean+init_mz_std*3];
y_coords_init = [y_lims_init(1), y_lims_init(1), y_lims_init(2), y_lims_init(2)];
x_patch1 = [time(1), time(moveidx(1)), time(moveidx(1)), time(1)];
patch(x_patch1, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
y_lims_mv = [mv_mz_mean-mv_mz_std*3, mv_mz_mean+mv_mz_std*3];
y_coords_init = [y_lims_mv(1), y_lims_mv(1), y_lims_mv(2), y_lims_mv(2)];
x_patch2 = [time(moveidx(2)), time(end), time(end), time(moveidx(2))];
patch(x_patch2, y_coords_init, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(time, rawData.mz);
plot(time(1:moveidx(1)), init_mz_mean*ones([moveidx(1), 1]), 'r--', 'LineWidth', 1);
plot(time(moveidx(2):end), mv_mz_mean*ones([length(time)-moveidx(2)+1, 1]), 'r--', 'LineWidth', 2);
y_lims = ylim;
y_coords= [y_lims(1), y_lims(1), y_lims(2), y_lims(2)];
x_patch = [time(moveidx(1)), time(moveidx(2)), time(moveidx(2)), time(moveidx(1))];
patch(x_patch, y_coords, 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('time (s)'); ylabel("mz");

init_mag = [init_mx_mean, init_my_mean, init_mz_mean, init_mx_std, init_my_std, init_mz_std];
mv_mag = [mv_mx_mean, mv_my_mean, mv_mz_mean, mv_mx_std, mv_my_std, mv_mz_std];
