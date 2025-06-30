close all; clear all; clc;
biasData = readtable("log\Huro_06_30_1432.csv");
    
bias_norm = NaN([length(biasData.debug0), 1]);
for i = 1:length(bias_norm)
    bias_norm(i) = norm([biasData.debug0(i), biasData.debug1(i), biasData.debug2(i)]);
end

ts = 1/570;
time = linspace(ts, ts*length(bias_norm), length(bias_norm));

rad2deg = 180/pi;
biasPlot = figure();
hold on; grid on;
% biasNorm = plot(time, bias_norm.*rad2deg, "LineWidth", 1);
% biasRol = plot(time, biasData.debug0.*rad2deg, "LineWidth", 1);
% biasPit = plot(time, biasData.debug1.*rad2deg, "LineWidth", 1);
% biasYaw = plot(time, biasData.debug2.*rad2deg, "LineWidth", 1);
biasNorm = plot(time, bias_norm , "LineWidth", 1);
biasRol = plot(time, biasData.debug0 , "LineWidth", 1);
biasPit = plot(time, biasData.debug1 , "LineWidth", 1);
biasYaw = plot(time, biasData.debug2 , "LineWidth", 1);
legend([biasNorm, biasRol, biasPit, biasYaw], {"Norm", "Rol", "Pit", "Yaw"});
