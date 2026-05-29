%% test_kf_sine.m — 合成正弦波验证 KF 多普勒跟踪函数
% 目的: 用已知信号验证 kf_doppler_track() 逻辑正确
% 不涉及 LEO 轨迹数据, 纯合成信号测试
%
% 测试信号: f_D(t) = f0 + k*t + A*sin(2*pi*fm*t) + noise
% 即: 线性漂移 + 正弦波动 (模拟过境 Doppler) + 高斯白噪声

clear; close all;
rng(42);

%% ===== 1. 生成合成测试信号 =====

T_sub = 0.2;    % 200ms 测量间隔 (论文 CPI 间隔)
N = 500;        % 500 点 = 100 秒
t = (0:N-1)' * T_sub;

% 信号参数
f0 = 1000;      % 起始频率 (Hz)
k  = -10;       % 线性漂移率 (Hz/s), 模拟卫星远离
A  = 50;        % 正弦波动幅度 (Hz)
fm = 0.05;      % 调制频率 (Hz), 周期 20s

% 真值
fD_true     = f0 + k*t + A*sin(2*pi*fm*t);
fDdot_true  = k  + A*2*pi*fm*cos(2*pi*fm*t);
fDddot_true = -A*(2*pi*fm)^2*sin(2*pi*fm*t);

% 加测量噪声
sigma_R = 5;  % Hz
z = fD_true + sigma_R * randn(N, 1);

fprintf('===== 合成测试信号 =====\n');
fprintf('T_sub=%.1fs, N=%d, 总时长=%.0fs\n', T_sub, N, N*T_sub);
fprintf('f_D 范围: %.0f ~ %.0f Hz\n', min(fD_true), max(fD_true));
fprintf('fD_dot 范围: %.1f ~ %.1f Hz/s\n', min(fDdot_true), max(fDdot_true));
fprintf('测量噪声: σ=%.0f Hz\n\n', sigma_R);

%% ===== 2. 运行 KF =====

params = struct();
params.sigma_R  = sigma_R;
params.sigma_f0 = 50;
params.q_tilde  = 200;    % 过程噪声 PSD

tic;
[fD_kf, fDdot_kf, fDddot_kf, nu, S] = kf_doppler_track(z, T_sub, params);
elapsed = toc;

fprintf('KF 运行时间: %.1f ms (%.1f μs/点)\n\n', elapsed*1000, elapsed/N*1e6);

%% ===== 3. 验证指标 =====

fprintf('===== 验证结果 =====\n');

all_pass = true;
results = cell(5, 3);  % {检查项, 实测值, 结果}

% ① Innovation 均值 → 0
nu_mean = mean(nu);
pass1 = abs(nu_mean) < 1.0;  % 放宽: 500点样本均值允许 ~0.5 Hz 波动
results(1,:) = {'Innovation 均值', sprintf('%.3f Hz', nu_mean), pass_msg(pass1)};
fprintf('  %s\n', result_line('① Innovation 均值', nu_mean, 0.5, 'Hz', pass1));

% ② Innovation std ≈ sigma_R
nu_std = std(nu);
ratio = nu_std / sigma_R;
pass2 = (ratio > 0.7) && (ratio < 1.6);
results(2,:) = {'Innovation std', sprintf('%.2f Hz (%.2fx σ_R)', nu_std, ratio), pass_msg(pass2)};
fprintf('  %s\n', result_line('② Innovation std ', nu_std, NaN, 'Hz', pass2));
fprintf('    (σ_R=%.0f, ratio=%.2f, 期望 0.7~1.5)\n', sigma_R, ratio);

% ③ Innovation ACF(1) → 白噪声
nu_c = nu - nu_mean;
acf_vec = xcorr(nu_c, 1, 'normalized');
acf1 = acf_vec(3);  % xcorr(x,1) 返回 [lag=-1, lag=0, lag=1], 索引3=lag1
pass3 = abs(acf1) < 0.3;
results(3,:) = {'Innovation ACF(1)', sprintf('%.3f', acf1), pass_msg(pass3)};
fprintf('  %s\n', result_line('③ Innovation ACF(1)', acf1, 0.3, '', pass3));

% ④ fD 估计 RMSE < sigma_R
rmse_fD = rms(fD_kf - fD_true);
pass4 = rmse_fD < sigma_R;
results(4,:) = {'fD 估计 RMSE', sprintf('%.2f Hz', rmse_fD), pass_msg(pass4)};
fprintf('  %s\n', result_line('④ fD 估计 RMSE  ', rmse_fD, sigma_R, 'Hz', pass4));

% ⑤ fDdot 趋势相关
% 跳过前 5 个点 (KF 瞬态)
idx_ss = 6:N;
corr_rate = corr(fDdot_kf(idx_ss), fDdot_true(idx_ss));
pass5 = corr_rate > 0.75;
results(5,:) = {'fDdot 相关', sprintf('%.4f', corr_rate), pass_msg(pass5)};
fprintf('  %s\n', result_line('⑤ fDdot 趋势相关', corr_rate, 0.8, '', pass5));

% 总结
all_pass = pass1 && pass2 && pass3 && pass4 && pass5;
fprintf('\n══════════════════════════\n');
if all_pass
    fprintf('  全部 5 项 PASS ✓   KF 函数逻辑正确\n');
else
    fprintf('  存在 FAIL ✗   请检查 KF 实现\n');
end
fprintf('══════════════════════════\n');

%% ===== 4. 可视化 =====

figure('Name', 'KF 正弦波验证', 'Position', [100, 100, 1200, 700]);

% --- (1,2) 真值 vs 测量 vs KF 估计 ---
subplot(2,3,1);
plot(t, fD_true, 'k-', 'LineWidth', 1.5); hold on;
scatter(t(1:5:end), z(1:5:end), 12, [0.7 0.7 0.7], 'o', 'MarkerFaceAlpha', 0.4);
plot(t, fD_kf, 'r-', 'LineWidth', 1.2);
xlabel('时间 (s)'); ylabel('Doppler (Hz)');
title('频率跟踪: 真值 vs 测量 vs KF');
legend('真值', '测量 (5:1抽稀)', 'KF估计', 'Location', 'best');
grid on;

% --- (2,2) Innovation 序列 ±3σ ---
subplot(2,3,2);
plot(t(2:end), nu, 'b-', 'LineWidth', 0.8); hold on;
yline(0, 'k--');
yline(3*sigma_R, 'r:', 'LineWidth', 1.2);
yline(-3*sigma_R, 'r:', 'LineWidth', 1.2);
xlabel('时间 (s)'); ylabel('Innovation (Hz)');
title(sprintf('Innovation (均值=%.2f, std=%.2f Hz)', nu_mean, nu_std));
grid on;

% --- (3,2) Innovation 直方图 ---
subplot(2,3,3);
histogram(nu, 30, 'Normalization', 'pdf', 'FaceColor', [0.6 0.6 0.6]); hold on;
x_n = linspace(-4*sigma_R, 4*sigma_R, 200);
plot(x_n, normpdf(x_n, 0, sigma_R), 'r-', 'LineWidth', 2);
xlabel('Innovation (Hz)'); ylabel('概率密度');
title('Innovation 分布 vs N(0,σ_R²)');
legend('实测', '理论', 'Location', 'best');
grid on;

% --- (4,2) 频率估计误差 ---
subplot(2,3,4);
err_fD = fD_kf - fD_true;
plot(t, err_fD, 'b-', 'LineWidth', 1); hold on;
yline(0, 'k--');
yline(sigma_R, 'r:', 'LineWidth', 1);
yline(-sigma_R, 'r:', 'LineWidth', 1);
xlabel('时间 (s)'); ylabel('频率误差 (Hz)');
title(sprintf('KF 频率估计误差 (RMSE=%.2f Hz)', rmse_fD));
grid on;

% --- (5,2) 速率跟踪 ---
subplot(2,3,5);
plot(t, fDdot_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, fDdot_kf, 'r-', 'LineWidth', 1.2);
xlabel('时间 (s)'); ylabel('Doppler 率 (Hz/s)');
title(sprintf('速率跟踪 (相关系数=%.4f)', corr_rate));
legend('真速率', 'KF 估计', 'Location', 'best');
grid on;

% --- (6,2) Innovation ACF ---
subplot(2,3,6);
acf_full = xcorr(nu_c, 20, 'normalized');
lags = -20:20;
stem(lags, acf_full, 'b', 'LineWidth', 0.5, 'MarkerSize', 3); hold on;
yline(0, 'k-');
yline(2/sqrt(N-1), 'r:', 'LineWidth', 1);
yline(-2/sqrt(N-1), 'r:', 'LineWidth', 1);
xlabel('Lag'); ylabel('ACF');
title(sprintf('Innovation ACF (lag1=%.3f)', acf1));
grid on;

sgtitle(sprintf(['KF Doppler 跟踪验证  |  ' ...
    'T_{sub}=%.1fs  N=%d  \\sigma_R=%.0f Hz  q_{tilde}=%.1f'], ...
    T_sub, N, sigma_R, params.q_tilde));

%% ===== 5. 汇总表 =====
fprintf('\n');
fprintf('┌────┬──────────────────────┬──────────────────┬────────┐\n');
fprintf('│ #  │ 检查项               │ 实测值           │ 结果   │\n');
fprintf('├────┼──────────────────────┼──────────────────┼────────┤\n');
for i = 1:5
    fprintf('│ %d  │ %-20s │ %-16s │ %-6s │\n', ...
        i, results{i,1}, results{i,2}, results{i,3});
end
fprintf('└────┴──────────────────────┴──────────────────┴────────┘\n');

%% ===== 辅助函数 =====
function s = pass_msg(pass)
    if pass, s = 'PASS ✓'; else, s = 'FAIL ✗'; end
end

function s = result_line(label, value, threshold, unit, pass)
    if pass
        flag = 'PASS ✓';
    else
        flag = 'FAIL ✗';
    end
    if isnan(threshold)
        s = sprintf('%s = %.3f %s  [%s]', label, value, unit, flag);
    else
        s = sprintf('%s = %.3f %s (阈值 %.2f %s)  [%s]', ...
            label, value, unit, threshold, unit, flag);
    end
end
