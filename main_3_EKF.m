%% main_3_EKF.m — EKF 实时定位 vs 多普勒噪声灵敏度
% 和 main_3.m 相同的数据 + 噪声扫描，但使用递推 EKF 代替批处理 NLS
% 验证: EKF 能否从相同噪声数据中逐步收敛到正确位置

clear all; close all;
rng(42);

%% ===== 参数配置 =====
f_c = 11.325e9;  c_light = 299792458;
el_threshold = 10;
max_iter_ekf = 3;  % IEKF 最大迭代次数 (处理非线性)

%% ===== 加载 0.2s 轨迹 =====
load('traj_3Sat_step02.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
T_sub = seconds(positionTT.Time(2) - positionTT.Time(1));
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

%% ===== 真实位置 & 真多普勒 =====
lat_rx = 31.82;  lon_rx = 117.23;  alt_rx = 30;
p_u_true = lla2ecef([lat_rx, lon_rx, alt_rx]);

doppler   = zeros(N_t, N_sat);
elevation = zeros(N_t, N_sat);

lat_r_rad = deg2rad(lat_rx);  lon_r_rad = deg2rad(lon_rx);
E_vec = [-sin(lon_r_rad), cos(lon_r_rad), 0];
N_vec = [-sin(lat_r_rad)*cos(lon_r_rad), -sin(lat_r_rad)*sin(lon_r_rad), cos(lat_r_rad)];
U_vec = [cos(lat_r_rad)*cos(lon_r_rad), cos(lat_r_rad)*sin(lon_r_rad), sin(lat_r_rad)];

for i = 1:N_sat
    for k = 1:N_t
        los  = satPos{i}(k,:) - p_u_true;
        dist = norm(los);
        v_r  = dot(satVel{i}(k,:), los) / dist;
        doppler(k,i) = -f_c / c_light * v_r;
        u = U_vec * los';  en = E_vec * los';  nn = N_vec * los';
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2 + nn^2)));
    end
end

fprintf('===== EKF 实时定位: 噪声灵敏度扫描 =====\n');
fprintf('轨迹: 3颗卫星, T_sub=%.1fs, N_t=%d\n', T_sub, N_t);
fprintf('Doppler 范围: %.0f ~ %.0f kHz\n\n', min(doppler(:))/1e3, max(doppler(:))/1e3);

%% ===== EKF 参数 =====
sigma_p0 = 50e3;      % 初始位置不确定度 50km
Q_diag   = (0.1)^2;   % 过程噪声 (m²), 静止接收机用极小值
P0 = diag([sigma_p0^2, sigma_p0^2, sigma_p0^2]);
Q  = diag([Q_diag, Q_diag, Q_diag]);

%% ===== 扫描多普勒噪声水平 =====
noise_levels = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500];  % Hz
N_levels = length(noise_levels);
N_trials = 50;

err_ekf_mean  = zeros(N_levels, 1);
err_ekf_std   = zeros(N_levels, 1);

fprintf('%-8s  %-14s  %-14s  %-10s  %-10s\n', ...
    'σ(Hz)', 'EKF误差(m)', 'NLS误差(m)', 'EKF收敛率', 'NLS收敛率');
fprintf('%-8s  %-14s  %-14s  %-10s  %-10s\n', ...
    '------', '----------', '----------', '--------', '--------');

for lev = 1:N_levels
    sigma_dop = noise_levels(lev);
    R = sigma_dop^2;
    trial_errs = zeros(N_trials, 1);
    n_conv = 0;

    for trial = 1:N_trials
        % --- 初始状态: 真值 + 50km 随机偏移 ---
        rng(trial);
        dp_init = sigma_p0 * randn(1, 3);
        x_post = p_u_true + dp_init;
        P_post = P0;

        % --- EKF 递推 ---
        for k = 1:N_t
            for sat_idx = 1:N_sat
                if elevation(k, sat_idx) < el_threshold
                    continue;
                end

                % 带噪声测量
                rng(lev * 1000 + trial * 100 + k * 10 + sat_idx);
                z = doppler(k, sat_idx) + sigma_dop * randn();

                p_s = satPos{sat_idx}(k, :)';
                v_s = satVel{sat_idx}(k, :)';

                % IEKF: 迭代重线性化
                for iter = 1:max_iter_ekf
                    l_vec = p_s - x_post;
                    dist  = norm(l_vec);
                    v_r   = dot(v_s, l_vec) / dist;
                    h_x   = -f_c / c_light * v_r;

                    % Jacobian (1×3), 复用 NLS 公式
                    v_dot_l = dot(v_s, l_vec);
                    l_dot_term = (v_dot_l / dist^2) * l_vec;
                    H = (f_c / (c_light * dist)) * (v_s - l_dot_term)';

                    % EKF Predict (首颗可见星时)
                    if iter == 1
                        x_pred = x_post;
                        P_pred = P_post + Q;
                    end

                    % EKF Update
                    S_k = H * P_pred * H' + R;
                    K   = P_pred * H' / S_k;
                    nu  = z - h_x;
                    x_post = x_pred + K * nu;
                    P_post = (eye(3) - K * H) * P_pred;
                end
            end
        end

        trial_errs(trial) = norm(x_post - p_u_true);
        if trial_errs(trial) < 500
            n_conv = n_conv + 1;
        end
    end

    err_ekf_mean(lev) = mean(trial_errs);
    err_ekf_std(lev)  = std(trial_errs);
    conv_rate_ekf = n_conv / N_trials * 100;

    % 对比 NLS 结果 (固定单次，同一数据)
    p_guess = p_u_true + sigma_p0 * randn(1, 3);
    err_nls = nls_position(p_guess, satPos, satVel, doppler, elevation, ...
        el_threshold, sigma_dop, f_c, c_light, p_u_true);

    fprintf('%-8.0f  %-14.2f  %-14.2f  %-9.0f%%  %-9.0f%%\n', ...
        sigma_dop, err_ekf_mean(lev), err_nls, conv_rate_ekf, 100);
end

%% ===== 单次详细演示 (σ=10Hz) =====
sigma_demo = 10;
R_demo = sigma_demo^2;

rng(1);
dp_init = sigma_p0 * randn(1, 3);
x_post = p_u_true + dp_init;
P_post = P0;

pos_history = zeros(N_t, 3);
err_history = zeros(N_t, 1);

meas_count = 0;
for k = 1:N_t
    for sat_idx = 1:N_sat
        if elevation(k, sat_idx) < el_threshold, continue; end

        rng(k * 10 + sat_idx);
        z = doppler(k, sat_idx) + sigma_demo * randn();

        p_s = satPos{sat_idx}(k, :)';
        v_s = satVel{sat_idx}(k, :)';

        for iter = 1:max_iter_ekf
            l_vec = p_s - x_post;
            dist  = norm(l_vec);
            h_x   = -f_c / c_light * dot(v_s, l_vec) / dist;
            v_dot_l = dot(v_s, l_vec);
            l_dot_term = (v_dot_l / dist^2) * l_vec;
            H = (f_c / (c_light * dist)) * (v_s - l_dot_term)';

            if iter == 1
                x_pred = x_post;
                P_pred = P_post + Q;
            end

            S_k = H * P_pred * H' + R_demo;
            K   = P_pred * H' / S_k;
            nu  = z - h_x;
            x_post = x_pred + K * nu;
            P_post = (eye(3) - K * H) * P_pred;
        end
        meas_count = meas_count + 1;
    end
    pos_history(k, :) = x_post';
    err_history(k) = norm(x_post - p_u_true);
end

fprintf('\n===== EKF 收敛演示 (σ=%.0f Hz) =====\n', sigma_demo);
fprintf('初始误差: %.1f km\n', norm(dp_init)/1e3);
fprintf('最终误差: %.1f m\n', err_history(end));
fprintf('总测量数: %d (%.0f 秒内, 仰角>%d°)\n', meas_count, N_t*T_sub, el_threshold);

%% ===== 可视化 =====

figure('Name', 'EKF 实时定位', 'Position', [100, 100, 1400, 500]);

% --- (1,3) 位置收敛过程 ---
subplot(1,3,1);
t_sec = (0:N_t-1)' * T_sub;
semilogy(t_sec, err_history, 'b-', 'LineWidth', 1.5);  hold on;
yline(10, 'r--', 'LineWidth', 1);
xlabel('时间 (s)');  ylabel('定位误差 (m)');
title(sprintf('EKF 收敛: %.0f km → %.1f m (σ=%.0f Hz)', ...
    norm(dp_init)/1e3, err_history(end), sigma_demo));
legend('3D 位置误差', '10m 线', 'Location', 'northeast');
grid on;

% --- (2,3) 定位轨迹 (ENU) ---
subplot(1,3,2);
% 转 ENU
pos_enu = zeros(N_t, 3);
true_enu = zeros(1, 3);
for k = 1:N_t
    pos_enu(k,:) = ecef2enu(pos_history(k,:), p_u_true, lat_rx, lon_rx, alt_rx);
end

% 每 50 点画一个标记看方向
step_plot = 50;
plot(pos_enu(1,1), pos_enu(1,2), 'rx', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(pos_enu(step_plot:step_plot:end,1), pos_enu(step_plot:step_plot:end,2), ...
    'b.-', 'MarkerSize', 8);
plot(0, 0, 'k+', 'MarkerSize', 15, 'LineWidth', 3);
xlabel('东向 (m)');  ylabel('北向 (m)');
title('EKF 位置估计轨迹 (ENU)');
legend('起始 (50km偏差)', 'EKF 收敛路径', '真值', 'Location', 'best');
grid on;  axis equal;

% --- (3,3) 噪声 vs 定位误差 ---
subplot(1,3,3);
errorbar(noise_levels, err_ekf_mean, err_ekf_std, 'b-o', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('多普勒噪声 σ (Hz)');  ylabel('定位误差 (m)');
title('EKF: 多普勒噪声 vs 定位误差');
grid on;

sgtitle(sprintf('EKF 实时定位 (T_{sub}=%.1fs, %d 颗卫星)', T_sub, N_sat));

%% ===== 辅助函数 =====

function err = nls_position(p_guess, satPos, satVel, doppler, elevation, ...
    el_threshold, sigma_dop, f_c, c_light, p_u_true)

    N_sat = length(satPos);
    N_t = size(doppler, 1);

    % 组装 NLS 输入
    P_all = [];  V_all = [];  y_all = [];
    for ii = 1:N_sat
        vis = elevation(:,ii) > el_threshold;
        P_all = [P_all; satPos{ii}(vis,:)];
        V_all = [V_all; satVel{ii}(vis,:)];
        y_all = [y_all; doppler(vis,ii)];
    end

    rng(1);
    y_noisy = y_all + sigma_dop * randn(size(y_all));

    p_u = p_guess;
    for iter = 1:50
        l_vec = P_all - p_u;
        dist = vecnorm(l_vec, 2, 2);
        v_r = sum(V_all .* l_vec, 2) ./ dist;
        f_model = -f_c / c_light * v_r;
        df = y_noisy - f_model;

        v_dot_l = sum(V_all .* l_vec, 2);
        l_dot_term = (v_dot_l ./ (dist.^2)) .* l_vec;
        H = (f_c ./ (c_light .* dist)) .* (V_all - l_dot_term);

        dp = H \ df;
        dp_norm = norm(dp);
        p_u = p_u + dp';
        if dp_norm < 1e-3, break; end
    end
    err = norm(p_u - p_u_true);
end

function enu = ecef2enu(xyz, ref_xyz, lat0, lon0, h0)
    % ECEF → ENU (East-North-Up)
    t = xyz(:)' - ref_xyz;
    lat_r = deg2rad(lat0);  lon_r = deg2rad(lon0);
    R = [-sin(lon_r),           cos(lon_r),           0;
         -sin(lat_r)*cos(lon_r), -sin(lat_r)*sin(lon_r), cos(lat_r);
          cos(lat_r)*cos(lon_r),  cos(lat_r)*sin(lon_r), sin(lat_r)];
    enu = (R * t')';
end
