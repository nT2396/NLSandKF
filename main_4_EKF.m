%% main_4_EKF.m — EKF 实时定位: 直线运动接收机 + LEO Doppler
% 接收机: 从合肥出发, 匀速向东 10m/s (36km/h)
% 卫星:  3 颗 Starlink, 0.2s 间隔轨迹
% EKF 状态: [x, y, z, vx, vy, vz]  6D 同时估计位置+速度
%
% 对比 main_3_EKF: 接收机从静止变成移动, 状态从 3D→6D

clear all; close all;
rng(42);

%% ===== 参数配置 =====
f_c = 11.325e9;  c_light = 299792458;
el_threshold = 10;
max_iter_ekf = 5;
sigma_dop = 10;   % Hz, 多普勒测量噪声

%% ===== 加载卫星轨迹 =====
load('traj_3Sat_step02.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
T_sub = seconds(positionTT.Time(2) - positionTT.Time(1));
t_sec = (0:N_t-1)' * T_sub;
selSat = [59, 74, 86];

satPos = cell(N_sat, 1);  satVel = cell(N_sat, 1);
for i = 1:N_sat
    posTT = positionTT.(sprintf('Sat%d', selSat(i)));
    velTT = velocityTT.(sprintf('Sat%d', selSat(i)));
    if iscell(posTT)
        pM = zeros(N_t, 3);  vM = zeros(N_t, 3);
        for k = 1:N_t, pM(k,:) = posTT{k};  vM(k,:) = velTT{k}; end
        satPos{i} = pM;  satVel{i} = vM;
    else
        satPos{i} = posTT;  satVel{i} = velTT;
    end
end

%% ===== 生成接收机直线运动轨迹 =====
lat_rx = 31.82;  lon_rx = 117.23;  alt_rx = 30;
p_start = lla2ecef([lat_rx, lon_rx, alt_rx]);    % 起点 (合肥)

% 东向单位向量 (ECEF)
e_east = [-sin(deg2rad(lon_rx)), cos(deg2rad(lon_rx)), 0];

v_true = 10;  % m/s, 向东匀速
p_u_true = zeros(N_t, 3);
v_u_true = zeros(N_t, 3);
for k = 1:N_t
    p_u_true(k, :) = p_start + v_true * t_sec(k) * e_east;
    v_u_true(k, :) = v_true * e_east;
end
dist_traveled = v_true * N_t * T_sub;
fprintf('===== 接收机直线运动 =====\n');
fprintf('起点: (%.3f, %.3f, %.3f) m (ECEF)\n', p_start);
fprintf('速度: %.0f m/s 向东,  总位移: %.1f km\n', v_true, dist_traveled/1e3);
fprintf('轨迹时长: %.0f s,  T_sub=%.1f s,  N_t=%d\n\n', N_t*T_sub, T_sub, N_t);

%% ===== 计算真多普勒 + 仰角 =====
doppler   = zeros(N_t, N_sat);
elevation = zeros(N_t, N_sat);

for i = 1:N_sat
    for k = 1:N_t
        los  = satPos{i}(k,:) - p_u_true(k,:);
        dist = norm(los);
        v_rel = satVel{i}(k,:) - v_u_true(k,:);
        v_r = dot(v_rel, los) / dist;
        doppler(k,i) = -f_c / c_light * v_r;

        % 仰角 (用起点算, 变化极小)
        u = dot([cos(deg2rad(lat_rx))*cos(deg2rad(lon_rx)), ...
                 cos(deg2rad(lat_rx))*sin(deg2rad(lon_rx)), ...
                 sin(deg2rad(lat_rx))], los);
        en = dot([-sin(deg2rad(lon_rx)), cos(deg2rad(lon_rx)), 0], los);
        nn = dot([-sin(deg2rad(lat_rx))*cos(deg2rad(lon_rx)), ...
                  -sin(deg2rad(lat_rx))*sin(deg2rad(lon_rx)), ...
                  cos(deg2rad(lat_rx))], los);
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2 + nn^2)));
    end
end

fprintf('Doppler 范围: %.0f ~ %.0f kHz\n\n', min(doppler(:))/1e3, max(doppler(:))/1e3);

%% ===== EKF 参数 =====
R = sigma_dop^2;       % 测量噪声协方差

% 初始不确定度: 模拟 GNSS 精度 (论文 Fig.6: GNSS可用期提供初始位置)
sigma_p0 = 10;         % GNSS 定位精度 ~10m
sigma_v0 = 1;          % 初始速度不确定度 ~1 m/s

P0 = blkdiag(sigma_p0^2 * eye(3), sigma_v0^2 * eye(3));

% 过程噪声: 离散白噪声加速度模型 (DWNA)
q_acc = 0.01;  % (m/s²)²/Hz, 加速度过程噪声 PSD (小值: 接收机基本匀速)
dt = T_sub;
Q_pos = [dt^3/3, 0, 0, dt^2/2, 0, 0;
          0, dt^3/3, 0, 0, dt^2/2, 0;
          0, 0, dt^3/3, 0, 0, dt^2/2;
          dt^2/2, 0, 0, dt, 0, 0;
          0, dt^2/2, 0, 0, dt, 0;
          0, 0, dt^2/2, 0, 0, dt];
Q = q_acc * Q_pos;

fprintf('===== EKF 参数 (论文 GNSS-aided 初始化) =====\n');
fprintf('sigma_p0 = %.0f m (GNSS 精度),  sigma_v0 = %.0f m/s\n', sigma_p0, sigma_v0);
fprintf('q_acc = %.2f (m/s²)²/Hz\n', q_acc);
fprintf('σ_Dop  = %.0f Hz\n\n', sigma_dop);

%% ===== EKF 递推 =====

% 初始状态: GNSS 最后定位 + 小偏差 (论文 Fig.6 GNSS可用期→初始状态)
rng(1);
x_post = [p_start + sigma_p0 * randn(1,3),  v_u_true(1,:) + sigma_v0 * randn(1,3)];  % 1×6
P_post = P0;

pos_est  = zeros(N_t, 3);
vel_est  = zeros(N_t, 3);
err_pos  = zeros(N_t, 1);
err_vel  = zeros(N_t, 1);
meas_cnt = 0;

for k = 1:N_t
    % ---- EKF Predict (每个时间步一次) ----
    F = [1, 0, 0, dt, 0,  0;
         0, 1, 0, 0,  dt, 0;
         0, 0, 1, 0,  0,  dt;
         0, 0, 0, 1,  0,  0;
         0, 0, 0, 0,  1,  0;
         0, 0, 0, 0,  0,  1];
    x_pred = (F * x_post')';   % 1×6
    P_pred = F * P_post * F' + Q;

    % ---- EKF Update: 每颗可见卫星 ----
    for sat_idx = 1:N_sat
        if elevation(k, sat_idx) < el_threshold, continue; end

        rng(k * 10 + sat_idx);
        z = doppler(k, sat_idx) + sigma_dop * randn();

        p_s = satPos{sat_idx}(k, :);   % 1×3
        v_s = satVel{sat_idx}(k, :);

        % IEKF 迭代
        x_upd = x_pred;
        for iter = 1:max_iter_ekf
            p_u = x_upd(1:3);   % 1×3
            v_u = x_upd(4:6);   % 1×3

            l_vec = p_s - p_u;          % 1×3
            dist  = norm(l_vec);
            v_rel = v_s - v_u;           % 1×3
            v_r   = dot(v_rel, l_vec) / dist;
            h_x   = -f_c / c_light * v_r;

            % Jacobian H (1×6)
            % ∂fD/∂p = fc/(c*d) * [v_rel - (v_rel·l)/d² * l]
            v_rel_dot_l = dot(v_rel, l_vec);
            l_term = (v_rel_dot_l / dist^2) * l_vec;   % 1×3
            H_p = (f_c / (c_light * dist)) * (v_rel - l_term);  % 1×3

            % ∂fD/∂v = fc/(c*d) * l  (接收机运动方向对 Doppler 的影响)
            H_v = (f_c / (c_light * dist)) * l_vec;    % 1×3

            H = [H_p, H_v];  % 1×6

            % Update
            S_k = H * P_pred * H' + R;
            K   = P_pred * H' / S_k;     % 6×1
            nu  = z - h_x;
            x_upd = x_pred + (K * nu)';   % 1×6
            P_pred = (eye(6) - K * H) * P_pred;
        end

        x_pred = x_upd;  % 为下一颗卫星更新预测值
        meas_cnt = meas_cnt + 1;
    end

    x_post = x_pred;
    P_post = P_pred;

    pos_est(k, :) = x_post(1:3);
    vel_est(k, :) = x_post(4:6);
    err_pos(k)    = norm(x_post(1:3) - p_u_true(k,:));
    err_vel(k)    = norm(x_post(4:6) - v_u_true(k,:));
end

fprintf('===== EKF 定位结果 (GNSS 初始化+LEO Doppler 维持) =====\n');
fprintf('初始位置误差: %.1f m (GNSS 精度)\n', err_pos(1));
fprintf('最终位置误差: %.1f m\n', err_pos(end));
fprintf('最终速度误差: %.2f m/s\n', err_vel(end));
fprintf('总测量次数:   %d (仰角>%d°)\n', meas_cnt, el_threshold);

%% ===== 可视化 =====

figure('Name', 'EKF 直线运动接收机定位', 'Position', [50, 50, 1400, 800]);

% --- (2,2) 位置收敛 ---
subplot(2,2,1);
semilogy(t_sec, err_pos, 'b-', 'LineWidth', 1.5);  hold on;
yline(10, 'r--', 'LineWidth', 1);
yline(100, 'r:', 'LineWidth', 0.8);
xlabel('时间 (s)');  ylabel('3D 位置误差 (m)');
title(sprintf('位置误差: %.1f m → %.1f m (GNSS→LEO维持)', err_pos(1), err_pos(end)));
legend('位置误差', '10m', '100m', 'Location', 'northeast');
grid on;

% --- (1,2) 速度收敛 ---
subplot(2,2,2);
plot(t_sec, vel_est(:,1), 'r-', 'LineWidth', 1);  hold on;
plot(t_sec, vel_est(:,2), 'g-', 'LineWidth', 1);
plot(t_sec, vel_est(:,3), 'b-', 'LineWidth', 1);
plot(t_sec, v_u_true(:,1), 'r--', 'LineWidth', 0.8);
plot(t_sec, v_u_true(:,2), 'g--', 'LineWidth', 0.8);
plot(t_sec, v_u_true(:,3), 'b--', 'LineWidth', 0.8);
xlabel('时间 (s)');  ylabel('速度 (m/s)');
title(sprintf('速度估计 (最终误差=%.2f m/s)', err_vel(end)));
legend('vx est', 'vy est', 'vz est', 'vx true', 'vy true', 'vz true');
grid on;

% --- (2,2) 位置估计 vs 真值 (ENU) ---
subplot(2,2,3);
% 转换到以起点为原点的 ENU
pos_enu_est = zeros(N_t, 2);
pos_enu_true = zeros(N_t, 2);
for k = 1:N_t
    dp_est = pos_est(k,:) - p_start;
    dp_true = p_u_true(k,:) - p_start;
    e_e = e_east;
    n_n = [-sin(deg2rad(lat_rx))*cos(deg2rad(lon_rx)), ...
           -sin(deg2rad(lat_rx))*sin(deg2rad(lon_rx)), ...
           cos(deg2rad(lat_rx))];
    pos_enu_est(k,1) = dot(dp_est, e_e);
    pos_enu_est(k,2) = dot(dp_est, n_n);
    pos_enu_true(k,1) = dot(dp_true, e_e);
    pos_enu_true(k,2) = dot(dp_true, n_n);
end

step_p = 50;
plot(pos_enu_est(1,1), pos_enu_est(1,2), 'rx', 'MarkerSize', 12, 'LineWidth', 2); hold on;
plot(pos_enu_est(step_p:step_p:end,1), pos_enu_est(step_p:step_p:end,2), ...
    'b.-', 'MarkerSize', 6);
plot(pos_enu_true(:,1), pos_enu_true(:,2), 'k-', 'LineWidth', 2);
plot(pos_enu_true(end,1), pos_enu_true(end,2), 'ks', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('东向 (m)');  ylabel('北向 (m)');
title('ENU 轨迹估计');
legend('EKF 起始 (GNSS)', 'EKF 轨迹', '真值(匀速向东)', '终点', 'Location', 'best');
grid on;  axis equal;

% --- (4,2) 3D 真实 vs 估计运动 ---
subplot(2,2,4);
% 卫星轨迹 + 接收机轨迹 (ECEF, 中心化)
p_center = mean(p_u_true, 1);
plot3((satPos{1}(:,1)-p_center(1))/1e3, ...
      (satPos{1}(:,2)-p_center(2))/1e3, ...
      (satPos{1}(:,3)-p_center(3))/1e3, 'Color', [0.7 0.7 0.7]); hold on;
plot3((satPos{2}(:,1)-p_center(1))/1e3, ...
      (satPos{2}(:,2)-p_center(2))/1e3, ...
      (satPos{2}(:,3)-p_center(3))/1e3, 'Color', [0.7 0.7 0.7]);
plot3((satPos{3}(:,1)-p_center(1))/1e3, ...
      (satPos{3}(:,2)-p_center(2))/1e3, ...
      (satPos{3}(:,3)-p_center(3))/1e3, 'Color', [0.7 0.7 0.7]);

plot3((p_u_true(:,1)-p_center(1))/1e3, ...
      (p_u_true(:,2)-p_center(2))/1e3, ...
      (p_u_true(:,3)-p_center(3))/1e3, 'k-', 'LineWidth', 2);
plot3((pos_est(end,1)-p_center(1))/1e3, ...
      (pos_est(end,2)-p_center(2))/1e3, ...
      (pos_est(end,3)-p_center(3))/1e3, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('X_{ECEF} (km)');  ylabel('Y_{ECEF} (km)');  zlabel('Z_{ECEF} (km)');
title('3D 视图: 卫星轨迹(灰) + 接收机轨迹(黑)');
legend('Sat59', 'Sat74', 'Sat86', '接收机真值', 'EKF 终点', 'Location', 'best');
grid on;  view(45, 30);

sgtitle(sprintf('EKF 实时定位: 直线运动接收机 (GNSS初始化, v=%.0f m/s, T_{sub}=%.1fs, σ_{Dop}=%.0f Hz)', ...
    v_true, T_sub, sigma_dop));
