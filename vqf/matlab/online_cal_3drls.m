close all; clear all; clc;
% rawData = readtable("log\cal_data\Huro_06_30_1545.csv");
rawData = readtable("log\cal_data\Huro_07_04_1048.csv");
input_format = 'MM.dd HH:mm:ss:SSS';
rawData.Time = datetime(rawData.Time, 'InputFormat', input_format);

raw_mx = rawData.mx;
raw_my = rawData.my;
raw_mz = rawData.mz;

data_len = length(raw_mx);

N = 9;      % cov matrix and input size
L = 2500;    % initial data length [L, N]

%%%%%%%%%%%%%%%%
init_L = zeros([N, N]);
init_R = zeros([N, 1]);
%%%%%%%%%%%%%%%%

count = 1;
init_flag = true;
hard_iron_log = NaN([data_len, 3]);
soft_iron_log = NaN([data_len, N]);
P_log = NaN([data_len, N*N]);
v_log = NaN([data_len, N]);
d_log = NaN([data_len, N]);
input_log = NaN([data_len, N]);
cal_mag_xyz = NaN([data_len, 3]);
est_param_log = NaN([data_len, N]);
output = 1;    % Fixed output
scale_factor = 16;
for i = 1:data_len
    mag_norm = norm([raw_mx(i), raw_my(i), raw_mz(i)]);
    raw_mx(i) = raw_mx(i)/scale_factor;
    raw_my(i) = raw_my(i)/scale_factor;
    raw_mz(i) = raw_mz(i)/scale_factor;
    
    input = input_vec([raw_mx(i), raw_my(i), raw_mz(i)]);
    input_log(i, :) = input;

    if count > L
        if init_flag
            est_param = init_L \ init_R;
            P = eye(N) / init_L;
            init_flag = false;
        end
        e_k = output - input * est_param;
        P = P - (P * input') / (1 + input * P * input') * (input * P);
        est_param = est_param + P * input' * e_k;
        est_param_log(i, :) = est_param;
    else
        init_R = init_R + input';
        init_L = init_L + input'*input;
    end
    if ~init_flag
        [hard, soft, v, d] = iron_log(est_param);
        cal_mag_xyz(i,:) = (soft*(v' * [raw_mx(i), raw_my(i), raw_mz(i)]' + hard'))';
        hard_iron_log(i,:) = hard;
        soft_iron_log(i,:) = reshape(soft', 1, []);
        v_log(i,:) = reshape(v', 1, []);
        d_log(i,:) = reshape(d', 1, []);
        P_log(i,:) = reshape(P', 1, []);
    end
    count = count + 1;
end
calPlot = figure();
subplot(1, 2, 1);
hold on; grid on;
scatter3(raw_mx, raw_my, raw_mz, 'Marker', '.');
axis equal;
subplot(1, 2, 2);
hold on; grid on;
scatter3(cal_mag_xyz(:,1), cal_mag_xyz(:,2), cal_mag_xyz(:,3), 'Marker', '.');
axis equal;

%% LSM results
mc_1 = raw_mx.^2;
mc_2 = raw_my.^2;
mc_3 = raw_mz.^2;
mc_4 = 2*raw_mx.*raw_my;
mc_5 = 2*raw_mx.*raw_mz;
mc_6 = 2*raw_my.*raw_mz;
mc_7 = raw_mx;
mc_8 = raw_my;
mc_9 = raw_mz;

% y : output 
y = ones(height(rawData), 1);

A = [mc_1 mc_2 mc_3 mc_4 mc_5 mc_6 mc_7 mc_8 mc_9];
x = (A' * A) \ (A' * y);
M = [x(1) x(4) x(5);
    x(4) x(2) x(6);
    x(5) x(6) x(3)];
b = [x(7), x(8), x(9)];

[V, D] = eig(M);

new_b = b * V;

% Scale the matrix to make it a unit sphere
scale_q = new_b(1)^2/(4*D(1, 1)) + new_b(2)^2/(4*D(2, 2)) + new_b(3)^2/(4*D(3, 3)) + 1;
scale_matrix = sqrt(D./scale_q);

hard_iron = [new_b(1)/D(1,1), new_b(2)/D(2, 2), new_b(3)/D(3, 3)]./2;
soft_iron = V*scale_matrix;
lsm_hard = reshape(hard_iron', 1, []);
lsm_soft = reshape(soft_iron', 1, []);
lsm_v = reshape(V', 1, []);

hardPlot = figure();
subplot(3, 1, 1);
hold on; grid on;
hard1 = plot(rawData.Time, hard_iron_log(:, 1));
lsm_hard1 = plot(rawData.Time, lsm_hard(1) * ones(data_len, 1), 'r--');
subplot(3, 1, 2);
hold on; grid on;
hard2 = plot(rawData.Time, hard_iron_log(:, 2));
lsm_hard2 = plot(rawData.Time, lsm_hard(2) * ones(data_len, 1),'r--');
subplot(3, 1, 3);
hold on; grid on;
hard3 = plot(rawData.Time, hard_iron_log(:, 3));
lsm_hard3 = plot(rawData.Time, lsm_hard(3) * ones(data_len, 1),'r--');

softPlot = figure();
subplot(3, 3, 1);
hold on; grid on;
soft1 = plot(rawData.Time, real(soft_iron_log(:, 1)));
lsm_soft1 = plot(rawData.Time, lsm_soft(1) * ones(data_len, 1),'r--');
subplot(3, 3, 2);
hold on; grid on;
soft2 = plot(rawData.Time, real(soft_iron_log(:, 2)));
lsm_soft2 = plot(rawData.Time, lsm_soft(2) * ones(data_len, 1),'r--');
subplot(3, 3, 3);
hold on; grid on;
soft3 = plot(rawData.Time, real(soft_iron_log(:, 3)));
lsm_soft3 = plot(rawData.Time, lsm_soft(3) * ones(data_len, 1),'r--');
subplot(3, 3, 4);
hold on; grid on;
soft4 = plot(rawData.Time, real(soft_iron_log(:, 4)));
lsm_soft4 = plot(rawData.Time, lsm_soft(4) * ones(data_len, 1),'r--');
subplot(3, 3, 5);
hold on; grid on;
soft5 = plot(rawData.Time, real(soft_iron_log(:, 5)));
lsm_soft5 = plot(rawData.Time, lsm_soft(5) * ones(data_len, 1),'r--');
subplot(3, 3, 6);
hold on; grid on;
soft6 = plot(rawData.Time, real(soft_iron_log(:, 6)));
lsm_soft6 = plot(rawData.Time, lsm_soft(6) * ones(data_len, 1),'r--');
subplot(3, 3, 7);
hold on; grid on;
soft7 = plot(rawData.Time, real(soft_iron_log(:, 7)));
lsm_soft7 = plot(rawData.Time, lsm_soft(7) * ones(data_len, 1),'r--');
subplot(3, 3, 8);
hold on; grid on;
soft8 = plot(rawData.Time, real(soft_iron_log(:, 8)));
lsm_soft8 = plot(rawData.Time, lsm_soft(8) * ones(data_len, 1),'r--');
subplot(3, 3, 9);
hold on; grid on;
soft9 = plot(rawData.Time, real(soft_iron_log(:, 9)));
lsm_soft9 = plot(rawData.Time, lsm_soft(9) * ones(data_len, 1),'r--');

vPlot = figure();
subplot(3, 3, 1);
hold on; grid on;
v1 = plot(rawData.Time, real(v_log(:, 1)));
lsm_v1 = plot(rawData.Time, lsm_v(1) * ones(data_len, 1), 'r--');
subplot(3, 3, 2);
hold on; grid on;
v2 = plot(rawData.Time, real(v_log(:, 2)));
lsm_v2 = plot(rawData.Time, lsm_v(2) * ones(data_len, 1), 'r--');
subplot(3, 3, 3);
hold on; grid on;
v3 = plot(rawData.Time, real(v_log(:, 3)));
lsm_v3 = plot(rawData.Time, lsm_v(3) * ones(data_len, 1), 'r--');
subplot(3, 3, 4);
hold on; grid on;
v4 = plot(rawData.Time, real(v_log(:, 4)));
lsm_v4 = plot(rawData.Time, lsm_v(4) * ones(data_len, 1), 'r--');
subplot(3, 3, 5);
hold on; grid on;
v5 = plot(rawData.Time, real(v_log(:, 5)));
lsm_v5 = plot(rawData.Time, lsm_v(5) * ones(data_len, 1), 'r--');
subplot(3, 3, 6);
hold on; grid on;
v6 = plot(rawData.Time, real(v_log(:, 6)));
lsm_v6 = plot(rawData.Time, lsm_v(6) * ones(data_len, 1), 'r--');
subplot(3, 3, 7);
hold on; grid on;
v7 = plot(rawData.Time, real(v_log(:, 7)));
lsm_v7 = plot(rawData.Time, lsm_v(7) * ones(data_len, 1), 'r--');
subplot(3, 3, 8);
hold on; grid on;
v8 = plot(rawData.Time, real(v_log(:, 8)));
lsm_v8 = plot(rawData.Time, lsm_v(8) * ones(data_len, 1), 'r--');
subplot(3, 3, 9);
hold on; grid on;
v9 = plot(rawData.Time, real(v_log(:, 9)));
lsm_v9 = plot(rawData.Time, lsm_v(9) * ones(data_len, 1), 'r--');

iron_logger = NaN([data_len, 1]);
for i = 1:length(iron_logger)
    iron_logger(i) = norm([soft_iron_log(i,:), hard_iron_log(i,:)]);
end

ironNormFig = figure();
hold on; grid on;
plot(rawData.Time, iron_logger);
ylabel("Norm of Soft/Hard iron");
%%
function input = input_vec(raw_mag)
    input = [raw_mag(1)^2, raw_mag(2)^2, raw_mag(3)^2, 2*raw_mag(1)*raw_mag(2), 2*raw_mag(1)*raw_mag(3), 2*raw_mag(2)*raw_mag(3), raw_mag(1), raw_mag(2), raw_mag(3)];
end

function [hard, soft, v, d] = iron_log(est_param)
    matrix = [est_param(1), est_param(4), est_param(5);...
              est_param(4), est_param(2), est_param(6);...
              est_param(5), est_param(6), est_param(3)];
    b = [est_param(7), est_param(8), est_param(9)];
    [V, D] = eig(matrix);
    new_b = b * V;
    scale_q = 1 + new_b(1)^2/(4*D(1, 1)) + new_b(2)^2/(4*D(2, 2)) + new_b(3)^2/(4*D(3, 3));
    scale_matrix = sqrt(D./scale_q);
    
    hard = [new_b(1)/D(1,1), new_b(2)/D(2, 2), new_b(3)/D(3, 3)]./2;
    soft = V*scale_matrix;
    v = V;
    d = D;
end
