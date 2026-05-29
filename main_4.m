%% main_4.m — 三阶 Kalman 滤波器: LEO 多普勒跟踪
% 论文: Shahcheraghi & Kassas, IEEE SPL 2024, Section III-C
%
% 状态向量: x = [f_D, fD_dot, fD_ddot]^T  (Hz, Hz/s, Hz/s²)
% 输入:     真多普勒 + 噪声 → 模拟 GLR/FFT 测量
% 输出:     KF 平滑多普勒 → NLS 定位
% 验证:     Innovation 白噪声检验 + 端到端定位精度

clear all; close all;
rng(42);

%% ==================== SECTION I: 加载数据 & 计算真值 ====================

% I.1 加载轨迹
load('traj_3Sat_step02.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
dt    = seconds(positionTT.Time(2) - positionTT.Time(1));
t_sec = (0:N_t-1)' * dt;

selSat = [59, 74, 86];
satNames = {'STARLINK-1270','STARLINK-1274','STARLINK-1309'}';

% 解析卫星位置速度
satPos = cell(N_sat,1); satVel = cell(N_sat,1);
for i = 1:N_sat
    posTT = positionTT.(sprintf('Sat%d', selSat(i)));
    velTT = velocityTT.(sprintf('Sat%d', selSat(i)));
    if iscell(posTT)
        posMat = zeros(N_t,3); velMat = zeros(N_t,3);
        for k = 1:N_t, posMat(k,:)=posTT{k}; velMat(k,:)=velTT{k}; end
        satPos{i}=posMat; satVel{i}=velMat;
    else
        satPos{i}=posTT; satVel{i}=velTT;
    end
end

% I.2 真实接收机位置 (合肥)
lat_rx=31.82; lon_rx=117.23; alt_rx=30;
p_u_true = lla2ecef([lat_rx, lon_rx, alt_rx]);

% I.3 计算真多普勒 + 多普勒率 + 多普勒加速度 + 仰角
f_c=11.325e9; c_light=299792458;
doppler    = zeros(N_t, N_sat);
doppler_rate = zeros(N_t, N_sat);
doppler_acc  = zeros(N_t, N_sat);
elevation  = zeros(N_t, N_sat);

lat_r_rad=deg2rad(lat_rx); lon_r_rad=deg2rad(lon_rx);
E_vec=[-sin(lon_r_rad), cos(lon_r_rad), 0];
N_vec=[-sin(lat_r_rad)*cos(lon_r_rad), -sin(lat_r_rad)*sin(lon_r_rad), cos(lat_r_rad)];
U_vec=[cos(lat_r_rad)*cos(lon_r_rad), cos(lat_r_rad)*sin(lon_r_rad), sin(lat_r_rad)];

for i = 1:N_sat
    for k = 1:N_t
        los = satPos{i}(k,:) - p_u_true;
        dist = norm(los);
        v_r = dot(satVel{i}(k,:), los) / dist;
        doppler(k,i) = -f_c / c_light * v_r;

        u=U_vec*los'; en=E_vec*los'; nn=N_vec*los';
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2+nn^2)));
    end
    % 数值微分求 rate 和 acc
    doppler_rate(:,i) = gradient(doppler(:,i), dt);
    doppler_acc(:,i)  = gradient(doppler_rate(:,i), dt);
end

% Jerk (Hz/s³): 用于调 q_tilde
doppler_jerk = gradient(doppler_acc, dt);
fprintf('实测 Doppler jerk RMS: %.2f Hz/s³\n', rms(doppler_jerk(:)));

%% ==================== SECTION II: KF 参数定义 ====================

% II.1 时间参数
T_sub = dt;  % 测量间隔 (s), 匹配轨迹采样

% II.2 状态转移矩阵 (连续 Wiener 过程加速度模型离散化)
F = [1, T_sub, T_sub^2/2;
     0, 1,     T_sub;
     0, 0,     1];

% II.3 过程噪声协方差
q_tilde = 1.0;  % Hz²/s⁵, 适配 T_sub=10s (rate 步间变化 ~26 Hz/s)
Q = q_tilde * [T_sub^5/20, T_sub^4/8,  T_sub^3/6;
               T_sub^4/8,  T_sub^3/3,  T_sub^2/2;
               T_sub^3/6,  T_sub^2/2,  T_sub];

% II.4 测量模型
H = [1, 0, 0];
sigma_R = 10;   % Hz, FFT 跟踪测量噪声
R = sigma_R^2;

% II.5 初始协方差 (反映 GLR 捕获精度)
sigma_f0    = 100;   % Hz,    GLR Doppler 精度
sigma_fd0   = 500;   % Hz/s,  rate 初始不确定度
sigma_fdd0  = 50;    % Hz/s², acc  初始不确定度
P0 = diag([sigma_f0^2, sigma_fd0^2, sigma_fdd0^2]);

fprintf('\nKF 参数:\n');
fprintf('  T_sub = %.0f s,  q_tilde = %.4f Hz²/s⁵\n', T_sub, q_tilde);
fprintf('  sigma_R = %.0f Hz,  sigma_f0 = %.0f Hz\n', sigma_R, sigma_f0);
fprintf('  Q(1,1)^0.5 (process f_D):   %.1f Hz\n', sqrt(Q(1,1)));
fprintf('  Q(2,2)^0.5 (process fD_dot): %.1f Hz/s\n', sqrt(Q(2,2)));
fprintf('  Q(3,3)^0.5 (process fD_ddot):%.2f Hz/s²\n', sqrt(Q(3,3)));

%% ==================== SECTION III: 逐卫星 KF 处理 ====================

kf_doppler  = zeros(N_t, N_sat);  % KF 平滑 Doppler
kf_rate     = zeros(N_t, N_sat);  % KF 估计 Doppler rate
kf_acc      = zeros(N_t, N_sat);  % KF 估计 Doppler acc
kf_nu       = zeros(N_t-1, N_sat);% Innovation
kf_S        = zeros(N_t-1, N_sat);% Innovation 协方差

for sat = 1:N_sat
    fD_true = doppler(:, sat);

    % III.1 初始状态: 模拟 GLR 捕获 (rate 从前两点差分估计)
    fD_dot_init = (fD_true(2) - fD_true(1)) / T_sub;
    x_post = [fD_true(1) + sigma_f0 * randn();
              fD_dot_init + sigma_fd0 * randn() / 5;
              0];
    P_post = P0;

    kf_doppler(1, sat) = x_post(1);
    kf_rate(1, sat)    = x_post(2);
    kf_acc(1, sat)     = x_post(3);

    % III.2 KF 循环
    for k = 1:N_t-1
        % --- Predict ---
        x_pred = F * x_post;
        P_pred = F * P_post * F' + Q;

        % --- 合成测量 (FFT 跟踪精度) ---
        z = fD_true(k+1) + sigma_R * randn();

        % --- Update ---
        nu = z - H * x_pred;
        S  = H * P_pred * H' + R;
        K  = P_pred * H' / S;
        x_post = x_pred + K * nu;
        P_post = (eye(3) - K * H) * P_pred;

        % --- 存储 ---
        kf_nu(k, sat) = nu;
        kf_S(k, sat)  = S;
        kf_doppler(k+1, sat) = x_post(1);
        kf_rate(k+1, sat)    = x_post(2);
        kf_acc(k+1, sat)     = x_post(3);
    end

    fprintf('  Sat %d (%s): KF 完成\n', selSat(sat), satNames{sat});
end

%% ==================== SECTION IV: Innovation 诊断 ====================

nu_all = kf_nu(:);
S_all  = kf_S(:);

% IV.1 基本统计
nu_mean = mean(nu_all);
nu_std  = std(nu_all);
fprintf('\n===== Innovation 诊断 =====\n');
fprintf('均值: %.2f Hz (期望 ≈ 0)\n', nu_mean);
fprintf('标准差: %.2f Hz (期望 ≈ %.0f Hz)\n', nu_std, sigma_R);

% IV.2 自相关 (lag=1)
nu_centered = nu_all - nu_mean;
[acf, lags] = xcorr(nu_centered, 1, 'normalized');
acf_lag1 = acf(lags == 1);
fprintf('ACF(lag=1): %.3f (期望 < 0.3)\n', acf_lag1);

% IV.3 NIS 检验
nis = nu_all.^2 ./ S_all;
nis_95 = mean(nis > chi2inv(0.95, 1));
fprintf('NIS 超 95%% 比例: %.1f%% (期望 ≈ 5%%)\n', nis_95*100);

% IV.4 Per-satellite Doppler 估计 RMSE
fprintf('\n===== Doppler 估计精度 =====\n');
fprintf('%-10s  %-12s  %-12s  %-12s\n', '卫星', '原始RMSE(Hz)', 'KF RMSE(Hz)', '改善比');
for sat = 1:N_sat
    raw_rmse = rms(sigma_R);  % 理论值
    kf_rmse  = rms(kf_doppler(:,sat) - doppler(:,sat));
 fprintf('Sat%d      %-12.1f  %-12.1f  %-12.1fx\n', ...
        selSat(sat), raw_rmse, kf_rmse, raw_rmse/kf_rmse);
end

%% ==================== SECTION V: NLS 端到端定位对比 ====================

fprintf('\n===== NLS 端到端定位对比 =====\n');

% V.1 组装 NLS 输入 (复用 main_2/3 逻辑)
el_threshold = 10;
P_sat_all = [];  V_sat_all = [];  y_true_all = [];

for ii = 1:N_sat
    vis_mask = elevation(:,ii) > el_threshold;
    P_sat_all = [P_sat_all; satPos{ii}(vis_mask,:)];
    V_sat_all = [V_sat_all; satVel{ii}(vis_mask,:)];
    y_true_all = [y_true_all; doppler(vis_mask,ii)];
end
NK_nls = size(P_sat_all, 1);

% V.2 合成 KF Doppler 向量 (与 P_sat 行序一致)
y_kf_all = zeros(NK_nls, 1);
row = 1;
for ii = 1:N_sat
    vis_mask = elevation(:,ii) > el_threshold;
    n_i = sum(vis_mask);
    y_kf_all(row:row+n_i-1) = kf_doppler(vis_mask, ii);
    row = row + n_i;
end

% V.3 对比: 多组噪声水平
noise_levels_nls = [10, 50, 100];
N_nls = length(noise_levels_nls);
err_raw_nls  = zeros(N_nls, 1);
err_kf_nls   = zeros(N_nls, 1);

for n = 1:N_nls
    sigma_n = noise_levels_nls(n);

    rng(n*100);
    y_noisy = y_true_all + sigma_n * randn(NK_nls, 1);
    err_raw_nls(n) = nls_solver(P_sat_all, V_sat_all, y_noisy, p_u_true);

    err_kf_nls(n) = nls_solver(P_sat_all, V_sat_all, y_kf_all, p_u_true);

    fprintf(' σ=%.0f Hz:  原始 → %.1f m,  KF → %.1f m  (%+.0f%%)\n', ...
        sigma_n, err_raw_nls(n), err_kf_nls(n), ...
        (err_raw_nls(n)-err_kf_nls(n))/err_raw_nls(n)*100);
end

%% ==================== SECTION VI: 可视化 ====================

figure('Name', '三阶 KF 多普勒跟踪', 'Position', [50, 50, 1400, 900]);
colors = lines(N_sat);

% --- (1,1) 真 Doppler vs 噪声测量 vs KF 估计 (Sat59 示例) ---
subplot(3,2,1);
sat_plot = 1;
plot(t_sec, doppler(:,sat_plot)/1e3, 'k-', 'LineWidth', 1.5); hold on;
rng(sat_plot);
z_plot = doppler(:,sat_plot) + sigma_R * randn(N_t,1);
scatter(t_sec, z_plot/1e3, 20, [0.7 0.7 0.7], 'o', 'MarkerFaceAlpha', 0.4);
plot(t_sec, kf_doppler(:,sat_plot)/1e3, 'r-', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('Doppler (kHz)');
title(sprintf('Sat%d (%s) — 真值 vs 测量 vs KF', selSat(sat_plot), satNames{sat_plot}));
legend('真值', '噪声测量 (\sigma=10Hz)', 'KF 估计', 'Location', 'best');
grid on;

% --- (1,2) Innovation 序列 ±3σ 界 ---
subplot(3,2,2);
for sat = 1:N_sat
    stairs(t_sec(2:end), kf_nu(:,sat), 'Color', colors(sat,:), 'LineWidth', 0.8); hold on;
end
yline(0, 'k--');
yline(3*sigma_R, 'r:', 'LineWidth', 1.5); yline(-3*sigma_R, 'r:', 'LineWidth', 1.5);
xlabel('时间 (s)'); ylabel('Innovation (Hz)');
title(sprintf('Innovation 序列 (±3σ=±%.0f Hz)', 3*sigma_R));
grid on;

% --- (2,1) Doppler 估计误差 ---
subplot(3,2,3);
for sat = 1:N_sat
    err_fD = kf_doppler(:,sat) - doppler(:,sat);
    plot(t_sec, err_fD, 'Color', colors(sat,:), 'LineWidth', 1.2); hold on;
end
yline(0, 'k--');
xlabel('时间 (s)'); ylabel('Doppler 误差 (Hz)');
title('KF 估计误差 vs 真值');
grid on;

% --- (2,2) 定位误差对比柱状图 ---
subplot(3,2,4);
bar_data = [err_raw_nls'; err_kf_nls']';
b = bar(noise_levels_nls, bar_data);
b(1).FaceColor = [0.6 0.6 0.6]; b(2).FaceColor = 'r';
xlabel('测量噪声 σ (Hz)'); ylabel('定位误差 (m)');
title('NLS 定位误差: 原始噪声 vs KF 平滑');
legend('原始噪声测量', 'KF 平滑 Doppler', 'Location', 'northwest');
grid on;
% 标注改善倍数
for n = 1:N_nls
    text(noise_levels_nls(n), max(err_raw_nls(n), err_kf_nls(n))*1.05, ...
        sprintf('%.1fx', err_raw_nls(n)/err_kf_nls(n)), ...
        'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
end

% --- (3,1) Innovation 直方图 ---
subplot(3,2,5);
histogram(nu_all, 25, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]); hold on;
x_norm = linspace(-4*sigma_R, 4*sigma_R, 200);
plot(x_norm, normpdf(x_norm, 0, sigma_R), 'r-', 'LineWidth', 2);
xlabel('Innovation (Hz)'); ylabel('概率密度');
title(sprintf('Innovation 分布 (均值=%.1f, std=%.1f)', nu_mean, nu_std));
legend('实测', sprintf('N(0,%.0f^2)', sigma_R));
grid on;

% --- (3,2) 三颗卫星 KF 跟踪总览 ---
subplot(3,2,6);
for sat = 1:N_sat
    plot(t_sec, doppler(:,sat)/1e3, '-', 'Color', colors(sat,:), 'LineWidth', 1.5); hold on;
    plot(t_sec, kf_doppler(:,sat)/1e3, '--', 'Color', colors(sat,:)*0.6, 'LineWidth', 1.2);
end
xlabel('时间 (s)'); ylabel('Doppler (kHz)');
title('三颗卫星 KF 跟踪总览 (实线:真值, 虚线:KF)');
legend(arrayfun(@(x) sprintf('Sat%d', selSat(x)), 1:N_sat, 'UniformOutput', false));
grid on;

sgtitle(sprintf(['三阶 KF 多普勒跟踪 (论文 Sec III-C)  |  ' ...
    'T_{sub}=%.0fs  q_{tilde}=%.4f  \\sigma_R=%.0f Hz'], T_sub, q_tilde, sigma_R));

%% ==================== 辅助函数: NLS 求解器 ====================

function err = nls_solver(P_sat, V_sat, y_meas, p_u_true)
    f_c = 11.325e9; c_light = 299792458;
    maxIter = 50; tol_dp = 1e-3;

    rng(1);
    delta_p = 50e3 * randn(1,3);
    p_u = p_u_true + delta_p;

    for iter = 1:maxIter
        l_vec  = P_sat - p_u;
        dist   = vecnorm(l_vec, 2, 2);
        v_r    = sum(V_sat .* l_vec, 2) ./ dist;
        f_model = -f_c / c_light * v_r;
        df = y_meas - f_model;

        v_dot_l = sum(V_sat .* l_vec, 2);
        l_dot_term = (v_dot_l ./ (dist.^2)) .* l_vec;
        H = (f_c ./ (c_light .* dist)) .* (V_sat - l_dot_term);

        dp = H \ df;
        dp_norm = norm(dp);
        p_u = p_u + dp';

        if dp_norm < tol_dp, break; end
    end
    err = norm(p_u - p_u_true);
end
