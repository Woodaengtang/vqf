clear all; close all; clc;
raw_data = readtable("log\Huro_06_24_1003.csv");
%%
raw_data_mx = raw_data.mx;
raw_data_my = raw_data.my;
raw_data_mz = raw_data.mz;

N = 9;
P = NaN([2*N, N]);
Y = NaN([N, 1]);
count = 1;
min_len = size(P, 1);
data_len = length(raw_data_mx);

loop_flag = true;
init_flag = true;
I = eye(N);
while loop_flag
    new_mag_x = raw_data_mx(count);
    new_mag_y = raw_data_my(count);
    new_mag_z = raw_data_mz(count);
    input = [new_mag_y^2, new_mag_z^2, new_mag_x*new_mag_y, new_mag_y*new_mag_z, new_mag_z*new_mag_x, new_mag_x, new_mag_y, new_mag_z, 1];
    output = -new_mag_x^2;
    if count > min_len*2
        if init_flag
            est_param = (P'*P)\P'*Y;
            P = I / (P'*P);
            init_flag = false;
        end
        e_k = output - input * est_param;
        P = P - (P * input') / (1 + input * P * input') * (input * P);
        est_param = est_param + P * input' * e_k;
    else
        P(count, :) = input;
        Y(count) = -new_mag_x.^2;
    end
    count = count + 1;

    if count > data_len
        loop_flag = false;
    end 
end

center = [         -2,   -est_param(3),   -est_param(5);...
         -est_param(3), -2*est_param(1),   -est_param(4);...
         -est_param(5),   -est_param(4), -2*est_param(2)]\[est_param(6); est_param(7); est_param(8)];

r = 35;
[sp_x, sp_y, sp_z] = sphere(100);
sp_x = sp_x*r;
sp_y = sp_y*r;
sp_z = sp_z*r;

rawMagXYZ = figure();
hold on; grid on;
% rawXYZ = scatter3(rawData.mx, rawData.my, rawData.mz, 'Marker', '.');
calXYZ = scatter3(raw_data_mx-center(1), raw_data_my-center(2), raw_data_mz-center(3), 'Marker', '.');
surf(sp_x, sp_y, sp_z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.1);
axis equal; 
xlabel('x'); ylabel('y'); zlabel('z');
view([-51, 32]);
