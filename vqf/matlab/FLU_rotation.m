clear; close all; clc;
origin = [0, 0, 0];
F_axis = [1, 0, 0];
L_axis = [0, 1, 0];
U_axis = [0, 0, 1];

duration = 3;
fps = 30;
total_frames = duration * fps;

angle_rpy = zeros([total_frames, 3]);
time = linspace(0, 5, total_frames/3);
angle_rpy(1:total_frames/3, 1) = 90*sin(2*pi*0.1*time)';
angle_rpy(total_frames/3+1:total_frames*2/3, 2) = 90*sin(2*pi*0.1*time)';
angle_rpy(total_frames*2/3+1:end, 3) = 90*sin(2*pi*0.1*time)';

fig1 = figure();
hold on;
grid on;
axis equal;
view(3);
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2]);
xlabel('X (World)');
ylabel('Y (World)');
zlabel('Z (World)');

q_F = quiver3(origin(1), origin(2), origin(3), F_axis(1), F_axis(2), F_axis(3), 1, 'r', 'LineWidth', 2);
q_L = quiver3(origin(1), origin(2), origin(3), L_axis(1), L_axis(2), L_axis(3), 1, 'g', 'LineWidth', 2);
q_U = quiver3(origin(1), origin(2), origin(3), U_axis(1), U_axis(2), U_axis(3), 1, 'b', 'LineWidth', 2);

t_F = text(F_axis(1), F_axis(2), F_axis(3), 'F', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold');
t_L = text(L_axis(1), L_axis(2), L_axis(3), 'L', 'Color', 'g', 'FontSize', 14, 'FontWeight', 'bold');
t_U = text(U_axis(1), U_axis(2), U_axis(3), 'U', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold');

h_title = title('FLU Coordinate System Rotation');
gif_name = 'asset/frame_rotation.gif';
for f = 1:total_frames
    % t = (frame - 1) / fps;
    % yaw   = 90 * sin(2 * pi * 0.1 * t);
    % % pitch = 20 * sin(2 * pi * 0.2 * t);
    % % roll  = 45 * sin(2 * pi * 0.5 * t);
    % pitch = 0;
    % roll  = 0;
    roll = angle_rpy(f, 1);
    pitch = angle_rpy(f, 2);
    yaw = angle_rpy(f, 3);

    eul = [deg2rad(yaw), deg2rad(pitch), deg2rad(roll)];
    rotm = eul2rotm(eul, 'ZYX');
    
    F_rotated = (rotm * F_axis')';
    L_rotated = (rotm * L_axis')';
    U_rotated = (rotm * U_axis')';
    
    set(q_F, 'UData', F_rotated(1), 'VData', F_rotated(2), 'WData', F_rotated(3));
    set(q_L, 'UData', L_rotated(1), 'VData', L_rotated(2), 'WData', L_rotated(3));
    set(q_U, 'UData', U_rotated(1), 'VData', U_rotated(2), 'WData', U_rotated(3));
    
    set(t_F, 'Position', F_rotated * 1.1);
    set(t_L, 'Position', L_rotated * 1.1);
    set(t_U, 'Position', U_rotated * 1.1);
    drawnow;
    pause(0.001);
    frame = getframe(fig1);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img,256);

    if f == 1
        imwrite(imind,cm,gif_name,'gif','Loopcount',1,'DelayTime',1/30);
    else
        imwrite(imind,cm,gif_name,'gif','WriteMode','append','DelayTime',1/30);
    end
end