%% 多普勒测量噪声 vs 定位误差 — 灵敏度分析
% 在不同多普勒估计误差下运行 NLS，绘制噪声-误差曲线
%
% 接收机: 合肥 (31.82°N, 117.23°E)
% 卫星:   Sat59, Sat74, Sat86

clear all;

%% ===== 参数配置 =====
f_c = 11.325e9;  c = 299792458;
el_threshold = 10;       % 仰角阈值 (度)
init_offset_km = 50;     % 初始猜测偏移 (km)
maxIter = 50;            % 最大迭代次数
tol_dp  = 1e-3;          % 收敛阈值: 位置修正量 < 1mm
N_trials = 100;           % 每种噪声水平的蒙特卡洛次数

%% ===== 加载星历 =====
load('traj_3Sat.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
selSat = [59, 74, 86];

satPos = cell(N_sat, 1);
satVel = cell(N_sat, 1);
for i = 1:N_sat
    posTT = positionTT.(sprintf('Sat%d', selSat(i)));
    velTT = velocityTT.(sprintf('Sat%d', selSat(i)));
    if iscell(posTT)
        posMat = zeros(N_t, 3); velMat = zeros(N_t, 3);
        for k = 1:N_t
            posMat(k,:) = posTT{k}; velMat(k,:) = velTT{k};
        end
        satPos{i} = posMat; satVel{i} = velMat;
    else
        satPos{i} = posTT; satVel{i} = velTT;
    end
end

%% ===== 真实位置 & 真实多普勒 =====
lat_rx = 31.82; lon_rx = 117.23; alt_rx = 30;
p_u_true = lla2ecef([lat_rx, lon_rx, alt_rx]);

doppler   = zeros(N_t, N_sat);
elevation = zeros(N_t, N_sat);
for i = 1:N_sat
    for k = 1:N_t
        los  = satPos{i}(k,:) - p_u_true;
        dist = norm(los);
        v_r  = dot(satVel{i}(k,:), los) / dist;
        doppler(k,i) = -f_c / c * v_r;

        lat_r_rad = deg2rad(lat_rx); lon_r_rad = deg2rad(lon_rx);
        los2 = satPos{i}(k,:) - p_u_true;
        E = [-sin(lon_r_rad), cos(lon_r_rad), 0];
        N = [-sin(lat_r_rad)*cos(lon_r_rad), -sin(lat_r_rad)*sin(lon_r_rad), cos(lat_r_rad)];
        U = [cos(lat_r_rad)*cos(lon_r_rad), cos(lat_r_rad)*sin(lon_r_rad), sin(lat_r_rad)];
        u = U * los2'; en = E * los2'; nn = N * los2';
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2 + nn^2)));
    end
end

%% ===== 组装 NLS 输入 (卫星位置/速度 + 无噪声多普勒) =====
P_sat = [];  V_sat = [];  y_meas_true = [];
for ii = 1:N_sat
    vis_mask = elevation(:,ii) > el_threshold;
    P_sat = [P_sat; satPos{ii}(vis_mask,:)];
    V_sat = [V_sat; satVel{ii}(vis_mask,:)];
    y_meas_true = [y_meas_true; doppler(vis_mask,ii)];
end
NK = size(P_sat, 1);
fprintf('总观测方程数: %d (卫星%d颗, 仰角>%d°)\n', NK, N_sat, el_threshold);

%% ===== 扫描不同多普勒噪声水平 =====
noise_levels = [0, 1, 2, 5, 10, 20, 50, 100, 200, 500];  % Hz (标准差)
N_levels = length(noise_levels);

% 存储结果
err_mean = zeros(N_levels, 1);   % 平均定位误差 (m)
err_std  = zeros(N_levels, 1);   % 定位误差标准差 (m)
err_all  = zeros(N_levels, N_trials);  % 所有 trial 结果

fprintf('\n===== 扫描多普勒噪声水平 =====\n');
fprintf('%-8s  %-12s  %-12s  %-8s\n', 'σ_Dop', '平均误差(m)', '误差std(m)', '收敛率');
fprintf('%-8s  %-12s  %-12s  %-8s\n', '------', '----------', '----------', '------');

for lev = 1:N_levels
    sigma_dop = noise_levels(lev);
    trial_errs = zeros(N_trials, 1);
    n_converged = 0;

    for trial = 1:N_trials
        % --- 加噪声 ---
        rng(lev * 1000 + trial);
        y_meas = y_meas_true + sigma_dop * randn(NK, 1);

        % --- 初始猜测 ---
        rng(trial);
        delta_p = init_offset_km * 1e3 * randn(1,3);
        p_u = p_u_true + delta_p;

        % --- NLS 迭代 ---
        for iter = 1:maxIter
            l_vec  = P_sat - p_u;
            dist   = vecnorm(l_vec, 2, 2);
            v_r    = sum(V_sat .* l_vec, 2) ./ dist;
            f_model = -f_c / c * v_r;
            df = y_meas - f_model;

            v_dot_l = sum(V_sat .* l_vec, 2);
            l_dot_term = (v_dot_l ./ (dist.^2)) .* l_vec;
            H = (f_c ./ (c .* dist)) .* (V_sat - l_dot_term);

            dp = H \ df;
            dp_norm = norm(dp);
            p_u = p_u + dp';

            if dp_norm < tol_dp
                n_converged = n_converged + 1;
                break;
            end
        end

        trial_errs(trial) = norm(p_u - p_u_true);
    end

    err_mean(lev) = mean(trial_errs);
    err_std(lev)  = std(trial_errs);
    conv_rate = n_converged / N_trials * 100;

    fprintf('%-8.0f  %-12.2f  %-12.2f  %-7.0f%%\n', ...
        sigma_dop, err_mean(lev), err_std(lev), conv_rate);
end

%% ===== 绘图 =====
figure('Name', '多普勒噪声 vs 定位误差', 'Position', [200, 200, 900, 400]);

% --- 子图1: 线性坐标 ---
subplot(1,2,1);
errorbar(noise_levels, err_mean, err_std, 'b-o', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b');
set(gca, 'XScale', 'linear', 'YScale', 'linear');
xlabel('多普勒估计误差 \sigma_{Dop} (Hz)');
ylabel('定位误差 (m)');
title('线性坐标');
grid on;

% 标注数值
for lev = 1:N_levels
    if err_mean(lev) < max(ylim)
        text(noise_levels(lev), err_mean(lev) + err_std(lev) + 2, ...
            sprintf('%.0f', err_mean(lev)), ...
            'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    end
end

% --- 子图2: log-log 坐标 ---
subplot(1,2,2);
idx_pos = noise_levels > 0;
x_log = noise_levels(idx_pos);
y_log = err_mean(idx_pos);
y_lo = max(err_mean(idx_pos) - err_std(idx_pos), 0.1);
y_hi = err_mean(idx_pos) + err_std(idx_pos);

loglog(x_log, y_log, 'r-s', ...
    'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r');
hold on;
% 误差棒 (逐段画)
for jj = 1:length(x_log)
    plot([x_log(jj), x_log(jj)], [y_lo(jj), y_hi(jj)], ...
        'r-', 'LineWidth', 0.5);
end
xlabel('多普勒估计误差 \sigma_{Dop} (Hz)');
ylabel('定位误差 (m)');
title('双对数坐标');
grid on;

sgtitle(sprintf('多普勒测量噪声 vs NLS定位误差 (初始偏差 %.0f km, %d trials)', ...
    init_offset_km, N_trials));

%% ===== 输出总结 =====
fprintf('\n===== 总结 =====\n');
fprintf('当多普勒估计误差 σ = 500 Hz 时:\n');
fprintf('  平均定位误差 = %.1f m\n', err_mean(end));
fprintf('  定位误差 std = %.1f m\n', err_std(end));
fprintf('当多普勒估计误差 σ = 10 Hz 时 (论文 KF 精度):\n');
fprintf('  平均定位误差 = %.1f m\n', err_mean(noise_levels == 10));
