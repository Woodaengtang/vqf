close all; clear all; clc;
readDiff = readtable("log\Diff_06_30_2331.csv");

input_format = 'MM.dd HH:mm:ss:SSS';
readDiff.Time = datetime(readDiff.Time, 'InputFormat', input_format);

for idx = 1 : length(readDiff.Time)
    if readDiff.Diff_Yaw(idx) < -180
        readDiff.Diff_Yaw(idx) = readDiff.Diff_Yaw(idx) + 360;
    end
    if readDiff.XSens_Yaw(idx) < 0
        readDiff.XSens_Yaw(idx) = readDiff.XSens_Yaw(idx) + 360;
    end
end

x_coord = [readDiff.Time(1), readDiff.Time(end)];
x_patch = [x_coord(1), x_coord(2), x_coord(2), x_coord(1)];
y_coord = [-2, 2];
y_patch = [y_coord(1), y_coord(1), y_coord(2), y_coord(2)];

plotDiff = figure();
plotDiff.Position = [680, 50, 914, 945];
subplot(3, 1, 1);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(readDiff.Time, readDiff.Diff_Roll);
xlabel("time(s)"); ylabel("\delta\phi (deg)"); subtitle("Diff between ICM w VQF6D and Xsens");
xlim(x_coord);
subplot(3, 1, 2);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(readDiff.Time, readDiff.Diff_Pitch);
xlabel("time(s)"); ylabel("\delta\theta (deg)");
xlim(x_coord);
subplot(3, 1, 3);
hold on; grid on;
patch(x_patch, y_patch, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(readDiff.Time, readDiff.Diff_Yaw);
xlabel("time(s)"); ylabel("\delta\psi (deg)");
xlim(x_coord);

plotAngle = figure();
plotAngle.Position = [680, 50, 914, 945];
subplot(3, 1, 1);
hold on; grid on;
plot(readDiff.Time, readDiff.Huro_Roll);
plot(readDiff.Time, readDiff.XSens_Roll);
xlabel("time(s)"); ylabel("\phi (deg)"); subtitle("ICM w VQF6D and Xsens");
xlim(x_coord);
subplot(3, 1, 2);
hold on; grid on;
plot(readDiff.Time, readDiff.Huro_Pitch);
plot(readDiff.Time, readDiff.XSens_Pitch);
xlabel("time(s)"); ylabel("\theta (deg)");
xlim(x_coord);
subplot(3, 1, 3);
hold on; grid on;
plot(readDiff.Time, readDiff.Huro_Yaw);
plot(readDiff.Time, readDiff.XSens_Yaw);
xlabel("time(s)"); ylabel("\psi (deg)");
xlim(x_coord);

