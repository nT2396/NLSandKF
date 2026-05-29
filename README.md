# NLSandKF — LEO 卫星多普勒定位：KF 跟踪 + NLS 定位

复现 Shahcheraghi & Kassas, *"A Computationally Efficient Approach for Acquisition and Doppler Tracking for PNT With LEO Megaconstellations,"* IEEE SPL, 2024.

## 完整链路

```
TLE 星历 → SGP4 轨道传播 → 几何多普勒(真值)
                               ↓
GLR 盲检测 (Automatic/)  →  初始 fD 估计
                               ↓
KF Doppler 跟踪 (3阶 Wiener 加速度模型) → 平滑 fD 序列
                               ↓
NLS 非线性最小二乘 → 接收机位置 (x, y, z)
```

## 文件说明

### 核心函数

| 文件 | 说明 |
|------|------|
| `kf_doppler_track.m` | **三阶 Kalman 滤波器**：状态 `[fD, fD_dot, fD_ddot]`，Wiener 过程加速度模型，可配置参数结构体 |
| `test_kf_sine.m` | KF 正弦波验证脚本：合成信号 5 项自动 PASS/FAIL 检查 |

### 主流程脚本（main_*.m）

| 文件 | 说明 | 输入 |
|------|------|------|
| `main.m` | 几何多普勒模型 + 仰角计算 | `traj_3Sat.mat` |
| `main_2.m` | NLS 单点定位验证（无噪声） | `traj_3Sat.mat` |
| `main_3.m` | 多普勒噪声 → 定位误差灵敏度扫描 | `traj_3Sat.mat` |
| `main_4.m` | KF 首次尝试（T_sub=10s, 失败记录） | `traj_3Sat.mat` |
| `main_5.m` | **Phase A 完整链路**：噪声 Doppler → KF 平滑 → NLS 定位 | `traj_3Sat_step02.mat` |

### 轨迹生成

| 文件 | 说明 |
|------|------|
| `gen_tle_v3.m` | 生成 3 颗卫星 0.2s 间隔轨迹 → `traj_3Sat_step02.mat` |
| `gen_tle_mat_v2.m` | 生成 4 颗卫星 10s 间隔轨迹 |
| `generate_tle_mat.m` | 生成全部卫星 10s 间隔轨迹 |

### 文档

| 文件 | 内容 |
|------|------|
| `NLS_流程梳理.md` | NLS 完整流程梳理（6 张流程图 + 数据维度流转） |
| `NLS_多普勒定位_方程推导.md` | NLS 方程推导（Taylor 线性化 → Jacobian → Gauss-Newton） |
| `Jacobian 的核心作用.md` | Jacobian 核心作用详解（物理含义 + 代码实现） |

### 数据文件

| 文件 | 说明 |
|------|------|
| `traj_3Sat.mat` | 3 颗卫星 10s 间隔轨迹（参考） |
| `traj_3Sat_step02.mat` | **3 颗卫星 0.2s 间隔轨迹**（main_5 使用） |
| `tle_starlink_2.txt` | 原始 TLE 两行轨道根数 |
| `tle_starlink.tle` | TLE 参考文件 |

## 关键设计决策

### 为什么 T_sub = 0.2s？

论文的 CPI（相干处理间隔）是 6.4ms，KF 子累积周期远小于 1s。最初使用 T_sub=10s 时 Innovation std 飙到 835 Hz（期望 ~10 Hz），根本原因是 LEO Doppler 率在 10s 内变化上千 Hz，Wiener 过程加速度模型的 Q 矩阵完全跟不上。0.2s 间隔下 Doppler 率变化仅几 Hz/s，KF 可以正常跟踪。

### q_tilde 数据驱动估计

Wiener 过程加速度模型：`d(fD_ddot)/dt = sqrt(q_tilde) · w(t)`。离散化后 `ΔfD_ddot ~ N(0, q_tilde·T_sub)`，因此 `q_tilde = var(ΔfD_ddot) / T_sub`。直接由轨迹 jerk 统计量计算，无需手动调参。

### KF 验证策略

先用**合成正弦波**（`fD = f0 + kt + Asin(2π·fm·t) + noise`）验证 KF 函数逻辑正确性（5 项检查），再用**真实轨迹数据**验证工程可用性。分离算法正确性和参数适配两个问题。

## 运行方式

```matlab
% 1. 生成轨迹（如需重新生成）
gen_tle_v3

% 2. 验证 KF 函数
test_kf_sine

% 3. 运行完整链路
main_5
```

## KF 函数接口

```matlab
[fD_kf, fDdot_kf, fDddot_kf, nu, S] = kf_doppler_track(z, T_sub, params);
```

| 参数 | 类型 | 默认值 | 含义 |
|------|------|--------|------|
| `params.sigma_R` | scalar | 10 | 测量噪声 std (Hz) |
| `params.sigma_f0` | scalar | 100 | 初始 Doppler 不确定度 (Hz) |
| `params.q_tilde` | scalar | 1.0 | 过程噪声 PSD (Hz²/s⁵) |
| `params.init_fD_dot` | scalar/[] | [] | 初始 Doppler 率 (空=自动差分) |

## 参考文献

1. Shahcheraghi & Kassas, "A Computationally Efficient Approach for Acquisition and Doppler Tracking for PNT With LEO Megaconstellations," IEEE SPL, 2024.
2. Neinavaie, Khalife & Kassas, "Acquisition, Doppler Tracking, and Positioning With Starlink LEO Satellites: First Results," IEEE TAES, 2022.
