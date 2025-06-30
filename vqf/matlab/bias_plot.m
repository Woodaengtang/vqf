close all; clear all; clc;
biasData = readtable("log\Huro_06_30_1222.csv");

gry_norm = NaN([length(biasData.debug0), 1]);
bias_norm = NaN([length(biasData.debug0), 1]);
for i = 1:length(bias_norm)
    bias_norm(i) = norm([biasData.debug0(i), biasData.debug1(i), biasData.debug2(i)]);
    gry_norm(i) = norm([biasData.gx(i), biasData.gy(i), biasData.gz(i)]);
end

ts = 1/570;
time = linspace(ts, ts*length(bias_norm), length(bias_norm));

rmse = [norm(mean(biasData.gx) - biasData.debug0), norm(mean(biasData.gy) - biasData.debug1), norm(mean(biasData.gz) - biasData.debug2)];

rad2deg = 180/pi;
biasPlot = figure();
biasPlot.Position = [1000, 106, 792, 1132];
subplot(4, 1, 1);
hold on; grid on;
normGryRaw = plot(time, gry_norm);
biasNorm = plot(time, bias_norm, 'LineWidth', 2, 'LineStyle', '--');
legend([normGryRaw, biasNorm], {"Raw data", "Estimated bias"});
ylabel("Gry norm");
subplot(4, 1, 2);
hold on; grid on;
gxRaw = plot(time, biasData.gx);
biasRol = plot(time, biasData.debug0, 'LineWidth', 2, 'LineStyle', '--');
ylabel("G_x");
subplot(4, 1, 3);
hold on; grid on;
gyRaw = plot(time, biasData.gy);
biasPit = plot(time, biasData.debug1, 'LineWidth', 2, 'LineStyle', '--');
ylabel("G_y");
subplot(4, 1, 4);
hold on; grid on;
gzRaw = plot(time, biasData.gz);
biasYaw = plot(time, biasData.debug2, 'LineWidth', 2, 'LineStyle', '--');
ylabel("G_z"); xlabel("time (s)");
