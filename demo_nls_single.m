%% Step 3: NLS 多普勒定位 — 框架
% 已知: 多普勒测量值 + 卫星星历 → 反解接收机位置
%
% 接收机: 合肥 (31.82°N, 117.23°E)
% 时间:   2026-05-25 08:38:00 ~ 08:42:00 UTC
% 卫星:   Sat59(STARLINK-1270), Sat74(STARLINK-1274), Sat86(STARLINK-1309)

clear all;

%% 3.1 加载星历
load('traj_3Sat.mat', 'positionTT', 'velocityTT');
N_sat = width(positionTT);
N_t   = height(positionTT);
dt    = seconds(positionTT.Time(2) - positionTT.Time(1));
t_sec = (0:N_t-1)' * dt;

% 解析卫星位置速度 
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

%% 3.2 设置真实接收机位置 (合肥)
lat_rx = 31.82;   % 北纬
lon_rx = 117.23;  % 东经
alt_rx = 30;      % 高度 (米, 合肥海拔约30m)
p_u_true = lla2ecef([lat_rx, lon_rx, alt_rx]);

%% 3.3 计算真实多普勒 (作为 NLS 的输入)
% 用真实接收机位置计算卫星到接收机的真实多普勒频率
% NLS 算法只知道这个结果, 不知道 p_u_true
f_c = 11.325e9;  c = 299792458;
doppler   = zeros(N_t, N_sat);
elevation = zeros(N_t, N_sat);

for i = 1:N_sat
    for k = 1:N_t
        los  = satPos{i}(k,:) - p_u_true;
        dist = norm(los);
        v_r  = dot(satVel{i}(k,:), los) / dist;
        doppler(k,i) = -f_c / c * v_r;
        % 计算仰角: ENU 坐标系下 u 分量
        lat_r = deg2rad(lat_rx); lon_r = deg2rad(lon_rx);
        los = satPos{i}(k,:) - p_u_true;
        E = [-sin(lon_r), cos(lon_r), 0];
        N = [-sin(lat_r)*cos(lon_r), -sin(lat_r)*sin(lon_r), cos(lat_r)];
        U = [cos(lat_r)*cos(lon_r), cos(lat_r)*sin(lon_r), sin(lat_r)];
        u = U * los';
        en = E * los'; nn = N * los';
        elevation(k,i) = rad2deg(atan2(u, sqrt(en^2 + nn^2)));
    end
end

%% 3.4 检查各卫星可见情况
N_use = N_sat;  % .mat 中只有4颗, 全部使用
fprintf('选取卫星: ');
for i = 1:N_use
    fprintf('Sat%d ', selSat(i));
end
fprintf('\n');

fprintf('\n各卫星可见时刻统计 (仰角>10°):\n');
for ii = 1:N_use
    vis_count = sum(elevation(:,ii) > 10);
    max_el = max(elevation(:,ii));
    fprintf('  Sat%d: 最大仰角=%.1f°, 可见时刻=%d/%d\n', selSat(ii), max_el, vis_count, N_t);
end

%% 3.5 组装 NLS 输入数据
% 把 N 颗卫星 × K 个时刻的数据"拉直"成向量
P_sat = [];  % NK × 3
V_sat = [];  % NK × 3
y_meas = []; % NK × 1

el_threshold = 10;  % 仰角阈值 (度)

for ii = 1:N_use
    vis_mask = elevation(:,ii) > el_threshold;
    P_sat = [P_sat; satPos{ii}(vis_mask,:)];     % K_i × 3
    V_sat = [V_sat; satVel{ii}(vis_mask,:)];     % K_i × 3
    y_meas = [y_meas; doppler(vis_mask,ii)];     % K_i × 1
end

% P 代表卫星位置 经过星历计算的    实际：未知
% V 代表速度                       实际：未知
% y_meas 代表多普勒频偏            实际：可由接收机计算出

NK = size(P_sat, 1);
fprintf('\n总观测方程数: %d (卫星%d颗)\n', NK, N_use);

%% 3.6 设置初始猜测位置
% 模拟"不知道接收机在哪"——在真实位置上加一个偏移
init_offset_km = 50;  % 初始误差约 50 km (3D向量模值)
rng(0);
delta_p = init_offset_km * 1e3 * randn(1,3);
p_u = p_u_true + delta_p;

fprintf('\n===== NLS 定位开始 =====\n');
fprintf('真实位置:   (%.2f°N, %.2f°E, %.0f m) → ECEF [%.1f, %.1f, %.1f] m\n', ...
    lat_rx, lon_rx, alt_rx, p_u_true);
lla_init = ecef2lla([p_u(1), p_u(2), p_u(3)]);
lat_init = lla_init(1); lon_init = lla_init(2);
fprintf('初始猜测:   (%.2f°N, %.2f°E) → ECEF [%.1f, %.1f, %.1f] m\n', ...
    lat_init, lon_init, p_u);
fprintf('初始误差:   %.1f km\n', norm(p_u - p_u_true)/1e3);

%% 3.7 NLS 迭代 (核心定位算法)
maxIter = 50;
tol     = 1e-9;  % 收敛阈值 (Hz)

fprintf('\n--- 迭代过程 ---\n');
fprintf('%-6s  %-14s  %-14s  %-14s\n', '迭代', '残差(Hz)', '位置修正(m)', '累计误差(km)');
fprintf('%-6s  %-14s  %-14s  %-14s\n', '----', '----------', '------------', '------------');

for iter = 1:maxIter
    % --- Step A: 用当前猜测位置算理论多普勒 ---
    l_vec = P_sat - p_u;            % NK × 3, 视线向量  卫星位置减去估计位置得到视线向量
    dist  = vecnorm(l_vec, 2, 2);   % NK × 1, 距离
    v_r   = sum(V_sat .* l_vec, 2) ./ dist;  % NK × 1, 径向速度
    f_model = -f_c / c * v_r;       % NK × 1, 在估计位置的理论多普勒，f_model代表有误差的估计
    
    % --- Step B: 计算残差 ---
    df = y_meas - f_model;  % NK × 1
    
    % --- Step C: 收敛判断 ---
    residual_norm = norm(df);
    cumul_err = norm(p_u - p_u_true) / 1e3;
    if iter == 1
        fprintf('%-6d  %-14.6e  %-14s  %-14.4f\n', iter, residual_norm, '---', cumul_err);
    end
    if residual_norm < tol
        fprintf('%-6d  %-14.6e  %-14.6e  %-14.4f  ← 收敛!\n', iter, residual_norm, 0, cumul_err);
        break;
    end
    
    % --- Step D: 计算 Jacobian ---
    % h = (f_c/(c*d)) * [v_s - (v_s·l/d²)·l]   (NK × 3)
    v_dot_l = sum(V_sat .* l_vec, 2);  % NK × 1, v_s · l
    l_dot_term = (v_dot_l ./ (dist.^2)) .* l_vec;  % NK × 3, (v_s·l/d²) * l  
    H = (f_c ./ (c .* dist)) .* (V_sat - l_dot_term); % NK × 3
    
    % --- Step E: 求解修正量 ---
    dp = H \ df;  % 3 × 1
    
    % --- Step F: 更新位置 ---
    dp_norm = norm(dp);
    fprintf('%-6d  %-14.6e  %-14.6e  %-14.4f\n', iter, residual_norm, dp_norm, cumul_err);
    p_u = p_u + dp';
end

if iter == maxIter
    fprintf('警告: 达到最大迭代次数 %d, 未收敛\n', maxIter);
end

%% 3.8 输出定位结果
p_u_est = p_u;
err = norm(p_u_est - p_u_true);
lla_est = ecef2lla([p_u_est(1), p_u_est(2), p_u_est(3)]);
lat_est = lla_est(1); lon_est = lla_est(2); alt_est = lla_est(3);

fprintf('\n===== 定位结果 =====\n');
fprintf('真实位置:   (%.4f°N, %.4f°E, %.1f m)\n', lat_rx, lon_rx, alt_rx);
fprintf('估计位置:   (%.4f°N, %.4f°E, %.1f m)\n', lat_est, lon_est, alt_est);
fprintf('定位误差:   %.2f m (%.4f km)\n', err, err/1e3);
fprintf('迭代次数:   %d\n', iter);
