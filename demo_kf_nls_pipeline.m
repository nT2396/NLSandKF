%% main_5.m — Phase A: KF+NLS 端到端仿真 (0.2s 轨迹)
% 链路: 真Doppler + 噪声 → KF 平滑 → NLS 定位
% 对比: 原始噪声 Doppler vs KF 平滑 Doppler 的定位精度
%
% 前置: traj_3Sat_step02.mat (gen_tle_v3.m 生成, sampleTime=0.2s)

clear all; close all;
rng(42);

%% ==================== SECTION I: 加载数据 & 计算真值 ====================

% I.1 加载 0.2s 间隔轨迹
load('traj_3Sat_step02.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
T_sub = seconds(positionTT.Time(2) - positionTT.Time(1));
t_sec = (0:N_t-1)' * T_sub;

selSat   = [59, 74, 86];
satNames = {'STARLINK-1270','STARLINK-1274','STARLINK-1309'}';

fprintf('===== 轨迹数据 =====\n');
fprintf('T_sub = %.2f s,  N_t = %d,  总时长 = %.0f s\n', T_sub, N_t, N_t*T_sub);

% I.2 解析卫星位置/速度
satPos = cell(N_sat, 1);  satVel = cell(N_sat, 1);
for i = 1:N_sat
    posTT = positionTT.(sprintf('Sat%d', selSat(i)));
    velTT = velocityTT.(sprintf('Sat%d', selSat(i)));
    if iscell(posTT)
        posMat = zeros(N_t, 3);  velMat = zeros(N_t, 3);
        for k = 1:N_t, posMat(k,:) = posTT{k};  velMat(k,:) = velTT{k}; end
        satPos{i} = posMat;  satVel{i} = velMat;
    else
        satPos{i} = posTT;  satVel{i} = velTT;
    end
end

% I.3 真实接收机位置 (合肥)
lat_rx = 31.82;  lon_rx = 117.23;  alt_rx = 30;
p_u_true = lla2ecef([lat_rx, lon_rx, alt_rx]);

% I.4 计算真多普勒 + 仰角
f_c = 11.325e9;  c_light = 299792458;

doppler      = zeros(N_t, N_sat);
doppler_rate = zeros(N_t, N_sat);
doppler_acc  = zeros(N_t, N_sat);
elevation    = zeros(N_t, N_sat);

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

        u  = U_vec * los';  en = E_vec * los';  nn = N_vec * los';
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2 + nn^2)));
    end
    doppler_rate(:,i) = gradient(doppler(:,i), T_sub);
    doppler_acc(:,i)  = gradient(doppler_rate(:,i), T_sub);
end

doppler_jerk = gradient(doppler_acc, T_sub);
fprintf('Doppler 范围: %.0f ~ %.0f kHz\n', min(doppler(:))/1e3, max(doppler(:))/1e3);
fprintf('Doppler rate RMS: %.1f Hz/s\n', rms(doppler_rate(:)));
fprintf('Doppler jerk  RMS: %.1f Hz/s^3\n', rms(doppler_jerk(:)));

%% ==================== SECTION II: KF 参数 ====================

% II.1 数据驱动 q_tilde: 从轨迹 jerk 估计
% Wiener 过程加速度模型: d(fD_ddot)/dt = sqrt(q_tilde) * w(t)
% ΔfD_ddot ~ N(0, q_tilde * T_sub), 所以 q_tilde = var(ΔfD_ddot) / T_sub
delta_acc = diff(doppler_acc, 1, 1);  % N_t-1 × N_sat
q_tilde_est = var(delta_acc(:)) / T_sub;
fprintf('\n===== KF 参数 =====\n');
fprintf('数据驱动 q_tilde = %.2f Hz²/s⁵ (来自轨迹 jerk RMS=%.1f Hz/s³)\n', ...
    q_tilde_est, rms(doppler_jerk(:)));

% II.2 测量噪声 (对应 GLR 输出精度)
sigma_R = 10;  % Hz

% II.3 KF 参数结构体
params = struct();
params.sigma_R  = sigma_R;
params.sigma_f0 = 100;   % GLR 初始 Doppler 不确定度
params.sigma_fd0 = 500;  % 初始 Doppler 率不确定度
params.sigma_fdd0 = 50;  % 初始加速度不确定度
params.q_tilde   = q_tilde_est;

% 显示过程噪声等效 std
Q = params.q_tilde * [T_sub^5/20, T_sub^4/8,  T_sub^3/6;
                       T_sub^4/8,  T_sub^3/3,  T_sub^2/2;
                       T_sub^3/6,  T_sub^2/2,  T_sub];
fprintf('σ_R = %.0f Hz\n', sigma_R);
fprintf('q_tilde = %.1f Hz²/s⁵\n', params.q_tilde);
fprintf('Q 过程噪声:  σ_fD=%.2f Hz,  σ_rate=%.2f Hz/s,  σ_acc=%.3f Hz/s²\n', ...
    sqrt(Q(1,1)), sqrt(Q(2,2)), sqrt(Q(3,3)));

%% ==================== SECTION III: KF 跟踪 (逐卫星) ====================

kf_doppler  = zeros(N_t, N_sat);
kf_rate     = zeros(N_t, N_sat);
kf_acc      = zeros(N_t, N_sat);
kf_nu       = zeros(N_t-1, N_sat);
kf_S        = zeros(N_t-1, N_sat);

for sat = 1:N_sat
    fD_true = doppler(:, sat);

    % 模拟 GLR 测量: 真 Doppler + 噪声
    rng(sat * 100);
    z = fD_true + sigma_R * randn(N_t, 1);

    % 初始 Doppler 率: 前两点差分 + 不确定度
    fD_dot_init = (fD_true(2) - fD_true(1)) / T_sub;
    params.init_fD_dot  = fD_dot_init + params.sigma_fd0 * randn() / 5;
    params.init_fD_ddot = 0;

    [fD_kf, fDdot_kf, fDddot_kf, nu, S] = kf_doppler_track(z, T_sub, params);

    kf_doppler(:, sat) = fD_kf;
    kf_rate(:, sat)    = fDdot_kf;
    kf_acc(:, sat)     = fDddot_kf;
    kf_nu(:, sat)      = nu;
    kf_S(:, sat)       = S;

    rmse_kf = rms(fD_kf - fD_true);
    fprintf('  Sat%d (%s): KF RMSE = %.2f Hz\n', selSat(sat), satNames{sat}, rmse_kf);
end

%% ==================== SECTION IV: Innovation 诊断 ====================

nu_all = kf_nu(:);
S_all  = kf_S(:);

nu_mean = mean(nu_all);
nu_std  = std(nu_all);

% ACF lag=1
nu_centered = nu_all - nu_mean;
[acf, lags] = xcorr(nu_centered, 1, 'normalized');
acf_lag1 = acf(lags == 1);

% NIS 检验
nis = nu_all.^2 ./ S_all;
nis_95 = mean(nis > chi2inv(0.95, 1));

fprintf('\n===== Innovation 诊断 =====\n');
fprintf('均值:     %+.2f Hz  (期望 ≈ 0)\n', nu_mean);
fprintf('标准差:   %.2f Hz  (期望 ≈ %.0f Hz)\n', nu_std, sigma_R);
fprintf('ACF(1):   %.3f     (期望 < 0.3)\n', acf_lag1);
fprintf('NIS>95%%:  %.1f%%    (期望 ≈ 5%%)\n', nis_95 * 100);

pass1 = abs(nu_mean) < 1.0;
pass2 = (nu_std / sigma_R > 0.7) && (nu_std / sigma_R < 1.6);
pass3 = abs(acf_lag1) < 0.3;
pass4 = nis_95 < 0.15;

fprintf('\n诊断结果: ');
if pass1 && pass2 && pass3 && pass4
    fprintf('全部 PASS ✓\n');
else
    fprintf('存在 FAIL ✗  (均值:%d std:%d ACF:%d NIS:%d)\n', pass1, pass2, pass3, pass4);
end

% Per-satellite RMSE
fprintf('\n===== Doppler 估计精度 =====\n');
fprintf('%-10s  %-14s  %-14s  %-10s\n', '卫星', '原始RMSE(Hz)', 'KF RMSE(Hz)', '改善比');
for sat = 1:N_sat
    raw_rmse = sigma_R;
    kf_rmse  = rms(kf_doppler(:,sat) - doppler(:,sat));
    fprintf('Sat%d      %-14.1f  %-14.1f  %-9.1fx\n', ...
        selSat(sat), raw_rmse, kf_rmse, raw_rmse/kf_rmse);
end

%% ==================== SECTION V: NLS 端到端定位 ====================

fprintf('\n===== NLS 端到端定位 =====\n');

% V.1 组装 NLS 输入 (仰角 > 10°)
el_threshold = 10;
P_sat_all = [];  V_sat_all = [];  y_true_all = [];
for ii = 1:N_sat
    vis_mask = elevation(:,ii) > el_threshold;
    P_sat_all   = [P_sat_all;   satPos{ii}(vis_mask,:)];
    V_sat_all   = [V_sat_all;   satVel{ii}(vis_mask,:)];
    y_true_all  = [y_true_all;  doppler(vis_mask,ii)];
end
NK_nls = size(P_sat_all, 1);
fprintf('NLS 观测方程数: %d (仰角 > %d°)\n', NK_nls, el_threshold);

% V.2 组装 KF Doppler (与 P_sat 行序一致)
y_kf_all = zeros(NK_nls, 1);
row = 1;
for ii = 1:N_sat
    vis_mask = elevation(:,ii) > el_threshold;
    n_i = sum(vis_mask);
    y_kf_all(row:row+n_i-1) = kf_doppler(vis_mask, ii);
    row = row + n_i;
end

% V.3 多噪声水平 Monte Carlo
noise_levels = [5, 10, 20, 50, 100];  % Hz
N_levels = length(noise_levels);
N_trials = 50;

err_noisy = zeros(N_levels, 1);  % 噪声 Doppler → NLS
err_kf    = zeros(N_levels, 1);  % KF 平滑 → NLS

fprintf('\n%-8s  %-14s  %-14s  %-10s\n', 'σ(Hz)', '噪声→NLS(m)', 'KF→NLS(m)', '改善比');
fprintf('%-8s  %-14s  %-14s  %-10s\n', '------', '----------', '----------', '------');

for n = 1:N_levels
    sigma_n = noise_levels(n);
    noisy_errs = zeros(N_trials, 1);

    for trial = 1:N_trials
        rng(n * 1000 + trial);
        y_noisy = y_true_all + sigma_n * randn(NK_nls, 1);
        noisy_errs(trial) = nls_solver(P_sat_all, V_sat_all, y_noisy, p_u_true);
    end
    err_noisy(n) = mean(noisy_errs);

    err_kf(n) = nls_solver(P_sat_all, V_sat_all, y_kf_all, p_u_true);

    fprintf('%-8.0f  %-14.1f  %-14.1f  %-9.1fx\n', ...
        sigma_n, err_noisy(n), err_kf(n), err_noisy(n)/err_kf(n));
end

%% ==================== SECTION VI: 可视化 ====================

figure('Name', 'KF+NLS 端到端仿真 (T_{sub}=0.2s)', 'Position', [50, 50, 1400, 900]);
colors = lines(N_sat);

% --- (1,2) 真值 vs 测量 vs KF (Sat59) ---
subplot(2,3,1);
sat_plot = 1;
plot(t_sec, doppler(:,sat_plot)/1e3, 'k-', 'LineWidth', 1.5); hold on;
rng(sat_plot*100);
z_plot = doppler(:,sat_plot) + sigma_R * randn(N_t,1);
scatter(t_sec(1:50:end), z_plot(1:50:end)/1e3, 8, [0.7 0.7 0.7], 'o', 'MarkerFaceAlpha', 0.3);
plot(t_sec, kf_doppler(:,sat_plot)/1e3, 'r-', 'LineWidth', 1.5);
xlabel('时间 (s)');  ylabel('Doppler (kHz)');
title(sprintf('Sat%d (%s) — 真值 vs 测量 vs KF', selSat(sat_plot), satNames{sat_plot}));
legend('真值', '噪声测量', 'KF 估计', 'Location', 'best');
grid on;

% --- (1,2) Innovation 序列 ±3σ ---
subplot(2,3,2);
for sat = 1:N_sat
    plot(t_sec(2:end), kf_nu(:,sat), 'Color', colors(sat,:), 'LineWidth', 0.6); hold on;
end
yline(0, 'k--');
yline(3*sigma_R, 'r:', 'LineWidth', 1.2);  yline(-3*sigma_R, 'r:', 'LineWidth', 1.2);
xlabel('时间 (s)');  ylabel('Innovation (Hz)');
title(sprintf('Innovation (均值=%.2f, std=%.1f, ACF(1)=%.3f)', nu_mean, nu_std, acf_lag1));
grid on;

% --- (2,1) Doppler 估计误差 ---
subplot(2,3,3);
for sat = 1:N_sat
    err_fD = kf_doppler(:,sat) - doppler(:,sat);
    plot(t_sec, err_fD, 'Color', colors(sat,:), 'LineWidth', 1); hold on;
end
yline(0, 'k--');  yline(sigma_R, 'r:', 'LineWidth', 1);  yline(-sigma_R, 'r:', 'LineWidth', 1);
xlabel('时间 (s)');  ylabel('Doppler 误差 (Hz)');
title('KF 估计误差');
legend(arrayfun(@(x) sprintf('Sat%d', selSat(x)), 1:N_sat, 'UniformOutput', false), ...
    'Location', 'best');
grid on;

% --- (2,2) NLS 定位误差对比 ---
subplot(2,3,4);
bar_data = [err_noisy, err_kf];
b = bar(noise_levels, bar_data);
b(1).FaceColor = [0.6 0.6 0.6];  b(2).FaceColor = 'r';
xlabel('测量噪声 σ (Hz)');  ylabel('定位误差 (m)');
title('NLS 定位: 噪声 vs KF 平滑');
legend('噪声 Doppler', 'KF 平滑', 'Location', 'northwest');
grid on;
for n = 1:N_levels
    text(noise_levels(n), max(err_noisy(n), err_kf(n)) * 1.08, ...
        sprintf('%.1fx', err_noisy(n)/err_kf(n)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
end

% --- (3,1) Innovation 直方图 ---
subplot(2,3,5);
histogram(nu_all, 40, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]); hold on;
x_norm = linspace(-4*sigma_R, 4*sigma_R, 200);
plot(x_norm, normpdf(x_norm, 0, sigma_R), 'r-', 'LineWidth', 2);
xlabel('Innovation (Hz)');  ylabel('概率密度');
title(sprintf('Innovation 分布 (std=%.1f Hz)', nu_std));
legend('实测', sprintf('N(0,%.0f²)', sigma_R), 'Location', 'best');
grid on;

% --- (3,2) 三颗卫星 KF 跟踪总览 ---
subplot(2,3,6);
for sat = 1:N_sat
    plot(t_sec, doppler(:,sat)/1e3, '-', 'Color', colors(sat,:), 'LineWidth', 1.2); hold on;
    plot(t_sec, kf_doppler(:,sat)/1e3, '--', 'Color', colors(sat,:)*0.6, 'LineWidth', 1);
end
xlabel('时间 (s)');  ylabel('Doppler (kHz)');
title('三颗卫星 KF 跟踪 (实线:真值, 虚线:KF)');
grid on;

sgtitle(sprintf('Phase A: KF+NLS 端到端 (T_{sub}=%.1fs  q_{tilde}=%.1f  \\sigma_R=%.0f Hz)', ...
    T_sub, params.q_tilde, sigma_R));

%% ==================== 辅助函数: NLS 求解器 ====================

function err = nls_solver(P_sat, V_sat, y_meas, p_u_true)
    f_c = 11.325e9;  c_light = 299792458;
    maxIter = 50;  tol_dp = 1e-3;

    % 初始猜测: 真值 + 50km 随机偏移
    rng(1);
    delta_p = 50e3 * randn(1, 3);
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
