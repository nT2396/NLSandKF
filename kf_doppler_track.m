function [fD_kf, fDdot_kf, fDddot_kf, nu, S] = kf_doppler_track(z, T_sub, params)
% kf_doppler_track  三阶 Kalman 滤波器: LEO 多普勒频率跟踪
%
% 论文: Shahcheraghi & Kassas, IEEE SPL 2024, Section III-C
% 状态: x = [f_D, fD_dot, fD_ddot]^T  (Hz, Hz/s, Hz/s²)
%
% 输入:
%   z       - N×1 多普勒测量序列 (Hz)
%   T_sub   - 测量间隔 (秒)
%   params  - KF 参数结构体 (可选字段, 未提供则用默认值)
%       .sigma_R    : 测量噪声 std (Hz),        默认 10
%       .sigma_f0   : 初始频率不确定度 (Hz),    默认 100
%       .sigma_fd0  : 初始速率不确定度 (Hz/s),  默认 500
%       .sigma_fdd0 : 初始加速度不确定度 (Hz/s²),默认 50
%       .q_tilde    : 过程噪声 PSD (Hz²/s⁵),    默认 1.0
%       .init_fD_dot: 手动指定初始速率 (Hz/s),  默认 [] (自动差分)
%       .init_fD_ddot:初始加速度 (Hz/s²),        默认 0
%
% 输出:
%   fD_kf     - N×1 KF 平滑 Doppler (Hz)
%   fDdot_kf  - N×1 KF 估计 Doppler 率 (Hz/s)
%   fDddot_kf - N×1 KF 估计 Doppler 加速度 (Hz/s²)
%   nu        - (N-1)×1 Innovation 序列 (Hz)
%   S         - (N-1)×1 Innovation 协方差 (Hz²)

    % ---- 默认参数 ----
    if nargin < 3, params = struct(); end
    sigma_R    = get_param(params, 'sigma_R',    10);
    sigma_f0   = get_param(params, 'sigma_f0',   100);
    sigma_fd0  = get_param(params, 'sigma_fd0',  500);
    sigma_fdd0 = get_param(params, 'sigma_fdd0', 50);
    q_tilde    = get_param(params, 'q_tilde',    1.0);
    init_fD_dot   = get_param(params, 'init_fD_dot',  []);
    init_fD_ddot  = get_param(params, 'init_fD_ddot', 0);

    N = length(z);
    if N < 3
        error('测量序列长度至少为 3');
    end

    % ---- 系统矩阵 ----
    F = [1, T_sub, T_sub^2/2;
         0, 1,     T_sub;
         0, 0,     1];

    Q = q_tilde * [T_sub^5/20, T_sub^4/8,  T_sub^3/6;
                   T_sub^4/8,  T_sub^3/3,  T_sub^2/2;
                   T_sub^3/6,  T_sub^2/2,  T_sub];

    H = [1, 0, 0];
    R = sigma_R^2;

    P0 = diag([sigma_f0^2, sigma_fd0^2, sigma_fdd0^2]);

    % ---- 初始化状态 ----
    if isempty(init_fD_dot)
        fD_dot_0 = (z(2) - z(1)) / T_sub;
    else
        fD_dot_0 = init_fD_dot;
    end

    x_post = [z(1); fD_dot_0; init_fD_ddot];
    P_post = P0;

    % ---- 预分配 ----
    fD_kf     = zeros(N, 1);
    fDdot_kf  = zeros(N, 1);
    fDddot_kf = zeros(N, 1);
    nu = zeros(N-1, 1);
    S  = zeros(N-1, 1);

    fD_kf(1)     = x_post(1);
    fDdot_kf(1)  = x_post(2);
    fDddot_kf(1) = x_post(3);

    % ---- KF 循环 ----
    for k = 1:N-1
        % Predict
        x_pred = F * x_post;
        P_pred = F * P_post * F' + Q;

        % Update
        nu_k   = z(k+1) - H * x_pred;
        S_k    = H * P_pred * H' + R;
        K_gain = P_pred * H' / S_k;
        x_post = x_pred + K_gain * nu_k;
        P_post = (eye(3) - K_gain * H) * P_pred;

        % 存储
        nu(k) = nu_k;
        S(k)  = S_k;
        fD_kf(k+1)     = x_post(1);
        fDdot_kf(k+1)  = x_post(2);
        fDddot_kf(k+1) = x_post(3);
    end
end

%% ===== 辅助函数 =====
function val = get_param(params, field, default)
    if isfield(params, field) && ~isempty(params.(field))
        val = params.(field);
    else
        val = default;
    end
end
