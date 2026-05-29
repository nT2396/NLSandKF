# NLS 多普勒定位 — 方程建立与求解推导

## 1. 问题描述

**已知**（可从测量/星历获得）：

| 符号 | 含义 | 来源 |
|------|------|------|
| $f_D^{(i)}(t_k)$ | 第 $i$ 颗卫星在 $t_k$ 时刻的多普勒频偏测量值 (Hz) | 接收机测量 |
| $\mathbf{p}_s^{(i)}(t_k) = [x_s, y_s, z_s]^T$ | 第 $i$ 颗卫星在 $t_k$ 时刻的 ECEF 位置 (m) | TLE + SGP4 |
| $\mathbf{v}_s^{(i)}(t_k) = [v_x, v_y, v_z]^T$ | 第 $i$ 颗卫星在 $t_k$ 时刻的 ECEF 速度 (m/s) | TLE + SGP4 |
| $f_c$ | 载波频率 (Hz)，Starlink = 11.325 GHz | 已知 |
| $c$ | 光速 = 299,792,458 m/s | 常数 |

**未知**：

| 符号 | 含义 | 维度 |
|------|------|------|
| $\mathbf{p}_u = [x_u, y_u, z_u]^T$ | 接收机 ECEF 位置 (m) | 3 |

---

## 2. 基础观测方程

对于卫星 $i$ 在时刻 $t_k$，多普勒频率与接收机位置的几何关系为：

$$\boxed{f_D^{(i)}(t_k) = -\frac{f_c}{c} \cdot \frac{\mathbf{v}_s^{(i)}(t_k) \cdot \big(\mathbf{p}_s^{(i)}(t_k) - \mathbf{p}_u\big)}{\big\|\mathbf{p}_s^{(i)}(t_k) - \mathbf{p}_u\big\|}}$$

其中：

- $\mathbf{l}^{(i)}(t_k) = \mathbf{p}_s^{(i)}(t_k) - \mathbf{p}_u$ 为视线向量 (Line-of-Sight)
- $d^{(i)}(t_k) = \|\mathbf{l}^{(i)}(t_k)\|$ 为星-地距离
- $v_r^{(i)}(t_k) = \dfrac{\mathbf{v}_s \cdot \mathbf{l}}{d}$ 为径向速度（接收机静止时 $\mathbf{v}_u = 0$）

**物理意义**：卫星靠近接收机时（$v_r < 0$），多普勒为正 → 频率升高；远离时（$v_r > 0$），多普勒为负 → 频率降低。

---

## 3. 方程数量分析

- **$N$ 颗卫星** × **$K$ 个时刻** → 共 $M = N \times K$ 个方程
- **未知数**：$\mathbf{p}_u = [x_u, y_u, z_u]^T$ → 3 个

| 场景 | 卫星数 | 时刻数 | 方程数 | 未知数 | 状态 |
|------|--------|--------|--------|--------|------|
| 单星单时刻 | 1 | 1 | 1 | 3 | 欠定（无穷多解） |
| 单星多时刻 | 1 | $K$ | $K$ | 3 | $K \ge 3$ 名义可解，但几何退化 |
| 多星单时刻 | $N$ | 1 | $N$ | 3 | $N \ge 3$ 名义可解，但瞬时几何差 |
| **多星多时刻** | **$N \ge 3$** | **$K \gg 1$** | **$NK$** | **3** | **超定 → NLS 求解** |

> **关键洞察**：多普勒是速率量，单时刻只能约束"接收机在某个等径向速度面上"，无法唯一确定位置。但卫星在运动，**多普勒曲线的形状**随接收机位置变化——不同位置看到的变化率（`df_D/dt`）不同，多时刻观测通过曲线形状差异来定位。

---

## 4. 非线性最小二乘 (NLS) 框架

### 4.1 代价函数

将 $M = NK$ 个观测方程组的残差平方和作为代价函数：

$$J(\mathbf{p}_u) = \sum_{i=1}^{N} \sum_{k=1}^{K} \Big[ f_D^{(i)}(t_k) - f_D^{\text{model}}(\mathbf{p}_u, t_k) \Big]^2 = \big\|\mathbf{f}_{\text{meas}} - \mathbf{f}_{\text{model}}(\mathbf{p}_u)\big\|^2$$

目标是找到 $\mathbf{p}_u$ 使 $J$ 最小：

$$\hat{\mathbf{p}}_u = \arg\min_{\mathbf{p}_u} J(\mathbf{p}_u)$$

### 4.2 Taylor 展开线性化

由于观测方程对 $\mathbf{p}_u$ 非线性（分母含距离 $d$），无法直接求解。在迭代点 $\hat{\mathbf{p}}_u^{(n)}$ 处做一阶 Taylor 展开：

$$\mathbf{f}_D(\mathbf{p}_u) \approx \mathbf{f}_D(\hat{\mathbf{p}}_u^{(n)}) + \underbrace{\frac{\partial \mathbf{f}_D}{\partial \mathbf{p}_u}\Big|_{\hat{\mathbf{p}}_u^{(n)}}}_{\mathbf{H}^{(n)}} \cdot (\mathbf{p}_u - \hat{\mathbf{p}}_u^{(n)})$$

记：
- $\Delta\mathbf{p} = \mathbf{p}_u - \hat{\mathbf{p}}_u^{(n)}$ — 位置修正量 (3×1)
- $\Delta\mathbf{f} = \mathbf{f}_{\text{meas}} - \mathbf{f}_D(\hat{\mathbf{p}}_u^{(n)})$ — 多普勒残差 ($M$×1)

得到**线性化方程组**：

$$\boxed{\mathbf{H}^{(n)} \cdot \Delta\mathbf{p} = \Delta\mathbf{f}}$$

其中 $\mathbf{H}^{(n)}$ 是 $M \times 3$ 的 Jacobian 矩阵。

---

## 5. Jacobian 矩阵推导

### 5.1 单行 Jacobian（一个卫星、一个时刻）

对于一个观测 $f_D^{(i)}(t_k)$，对接收机位置 $\hat{\mathbf{p}}_u$ 求偏导。

令 $\mathbf{l} = \mathbf{p}_s - \hat{\mathbf{p}}_u$（视线向量，1×3），$d = \|\mathbf{l}\|$（距离）。

**重要**：$\partial\mathbf{l}/\partial\hat{\mathbf{p}}_u = -\mathbf{I}_{3\times3}$（因为 $\mathbf{l} = \mathbf{p}_s - \hat{\mathbf{p}}_u$）

#### Step 1：径向速度的偏导

$$v_r = \frac{\mathbf{v}_s \cdot \mathbf{l}}{d}$$

使用商法则：

$$\frac{\partial v_r}{\partial\hat{\mathbf{p}}_u} = \frac{\frac{\partial}{\partial\hat{\mathbf{p}}_u}(\mathbf{v}_s \cdot \mathbf{l}) \cdot d - (\mathbf{v}_s \cdot \mathbf{l}) \cdot \frac{\partial d}{\partial\hat{\mathbf{p}}_u}}{d^2}$$

其中：
- $\frac{\partial}{\partial\hat{\mathbf{p}}_u}(\mathbf{v}_s \cdot \mathbf{l}) = \mathbf{v}_s \cdot \frac{\partial\mathbf{l}}{\partial\hat{\mathbf{p}}_u} = -\mathbf{v}_s^T$ (1×3)
- $\frac{\partial d}{\partial\hat{\mathbf{p}}_u} = \frac{\mathbf{l}^T}{d} \cdot \frac{\partial\mathbf{l}}{\partial\hat{\mathbf{p}}_u} = -\frac{\mathbf{l}^T}{d}$ (1×3)

代入：

$$\frac{\partial v_r}{\partial\hat{\mathbf{p}}_u} = \frac{(-\mathbf{v}_s^T) \cdot d - (\mathbf{v}_s \cdot \mathbf{l}) \cdot (-\frac{\mathbf{l}^T}{d})}{d^2} = -\frac{\mathbf{v}_s^T}{d} + \frac{(\mathbf{v}_s \cdot \mathbf{l})}{d^3} \cdot \mathbf{l}^T$$

#### Step 2：多普勒的偏导

$$f_D = -\frac{f_c}{c} \cdot v_r$$

$$\frac{\partial f_D}{\partial\hat{\mathbf{p}}_u} = -\frac{f_c}{c} \cdot \left(-\frac{\mathbf{v}_s^T}{d} + \frac{(\mathbf{v}_s \cdot \mathbf{l})}{d^3} \cdot \mathbf{l}^T\right)$$

$$\boxed{\frac{\partial f_D}{\partial\hat{\mathbf{p}}_u} = \frac{f_c}{c \cdot d} \left[ \mathbf{v}_s^T - \frac{\mathbf{v}_s \cdot \mathbf{l}}{d^2} \cdot \mathbf{l}^T \right]}$$

这就是 Jacobian 的**一行**（1×3 向量），对应 MATLAB 代码：

```matlab
v_dot_l    = dot(v_s, l);                        % 标量: v_s · l
l_dot_term = (v_dot_l / dist^2) * l;             % 1×3: (v_s·l/d²)·l
H_row      = (fc / (c * dist)) * (v_s - l_dot_term); % 1×3: Jacobian行
```

### 5.2 完整 Jacobian 矩阵

将所有 $M = NK$ 个观测堆叠：

$$\mathbf{H} = \begin{bmatrix}
\frac{\partial f_D^{(1)}(t_1)}{\partial x_u} & \frac{\partial f_D^{(1)}(t_1)}{\partial y_u} & \frac{\partial f_D^{(1)}(t_1)}{\partial z_u} \\[4pt]
\frac{\partial f_D^{(1)}(t_2)}{\partial x_u} & \frac{\partial f_D^{(1)}(t_2)}{\partial y_u} & \frac{\partial f_D^{(1)}(t_2)}{\partial z_u} \\[4pt]
\vdots & \vdots & \vdots \\[4pt]
\frac{\partial f_D^{(N)}(t_K)}{\partial x_u} & \frac{\partial f_D^{(N)}(t_K)}{\partial y_u} & \frac{\partial f_D^{(N)}(t_K)}{\partial z_u}
\end{bmatrix}_{M \times 3}$$

MATLAB 向量化实现（`main_2.m` 第144-148行）：

```matlab
l_vec = P_sat - p_u;                      % M×3: 所有视线向量 (M = NK)
dist  = vecnorm(l_vec, 2, 2);             % M×1: 所有距离
v_dot_l = sum(V_sat .* l_vec, 2);         % M×1: v_s · l (逐行点积)
l_dot_term = (v_dot_l ./ (dist.^2)) .* l_vec; % M×3: (v_s·l/d²)·l
H = (f_c ./ (c * dist)) .* (V_sat - l_dot_term); % M×3: 完整 Jacobian
```

---

## 6. 求解修正量

### 6.1 最小二乘解

$M \gg 3$ 时，$\mathbf{H}\Delta\mathbf{p} = \Delta\mathbf{f}$ 是超定方程组（方程数 > 未知数），用最小二乘求解：

$$\Delta\hat{\mathbf{p}} = (\mathbf{H}^T\mathbf{H})^{-1}\mathbf{H}^T \Delta\mathbf{f}$$

在 MATLAB 中，反斜杠运算符自动选择最优方法：

```matlab
dp = H \ df;   % 3×1: 最小二乘解 (MATLAB自动用QR分解)
```

等价于 `dp = (H'*H) \ (H'*df)` 或 `dp = inv(H'*H) * (H'*df)`。

### 6.2 加权的变体

如果不同卫星/不同时刻的测量精度不同（如低仰角噪声大），可引入权重矩阵 $\mathbf{W}$：

$$\Delta\hat{\mathbf{p}} = (\mathbf{H}^T\mathbf{W}\mathbf{H})^{-1}\mathbf{H}^T\mathbf{W} \Delta\mathbf{f}$$

其中 $\mathbf{W} = \text{diag}(w_1, w_2, \ldots, w_M)$，权重可设为与仰角相关（高仰角 = 高权重）。

---

## 7. Gauss-Newton 迭代

```
初始猜测 p_u^(0)
for n = 0, 1, 2, ... until convergence:
    1. l = P_sat - p_u^(n)                # 计算视线向量
    2. d = ||l||                          # 计算距离
    3. f_model = -(fc/c) * (V_sat·l)/d    # 理论多普勒
    4. df = f_meas - f_model              # 残差 (M×1)
    5. H = (fc/(c*d)) * [V_sat - (V_sat·l)/d² * l]  # Jacobian (M×3)
    6. dp = H \ df                        # 最小二乘修正量 (3×1)
    7. p_u^(n+1) = p_u^(n) + dp^T         # 更新位置
    8. if ||dp|| < ε: 停止                 # 收敛判断
```

**收敛判断**（论文条件）：位置修正量范数 $\|\Delta\mathbf{p}\| < \epsilon$（如 $10^{-6}$ m 即 1 微米），或者残差变化 $|J^{(n+1)} - J^{(n)}| < 10^{-9}$。

### 7.1 代码对应

| 步骤 | `main_2.m` 行号 |
|------|----------------|
| 1-3 | 第125-128行 (`l_vec`, `dist`, `v_r`, `f_model`) |
| 4 | 第131行 (`df`) |
| 5 | 第146-148行 (`H`) |
| 6 | 第151行 (`dp = H \ df`) |
| 7 | 第156行 (`p_u = p_u + dp'`) |
| 8 | 第139行 (当前用 `norm(df)`, 建议改为 `norm(dp)`) |

---

## 8. 数值示例

以合肥接收机 (31.82°N, 117.23°E) 为例：

| 参数 | 值 |
|------|-----|
| 卫星数 $N$ | 3 (Sat59, Sat74, Sat86) |
| 时间跨度 | 300 秒 (08:38–08:42 UTC) |
| 采样间隔 | 1 秒 |
| 总观测数 $M$ | 3 × 300 = 900 |
| 未知数 | 3 (x, y, z) |
| 超定比 $M/3$ | 300× |
| 初始偏差 | 50 km |
| 载频 $f_c$ | 11.325 GHz |

**单行 Jacobian 的物理量级**（Starlink, ~500km 轨道）：

| 分量 | 量级 |
|------|------|
| $d$ (距离) | ~500–2000 km |
| $\mathbf{v}_s$ (卫星速度) | ~7.5 km/s |
| $f_c/(c \cdot d)$ | ~0.075 Hz/m² |
| $\partial f_D/\partial x_u$ (单行) | ~0.1–1 Hz/m |

> 即接收机移动 1 米，多普勒变化约 0.1–1 Hz。Jacobian 矩阵的**条件数**取决于卫星几何分布——卫星散布越广，条件数越小，求解越稳定。

---

## 9. 扩展到含钟漂的模型

实际接收机的多普勒测量还包含卫星钟漂和接收机钟漂的影响。扩展后的模型：

### 9.1 含钟漂的观测方程

$$f_D^{(i)}(t_k) = -\frac{f_c}{c} \cdot \frac{\mathbf{v}_s^{(i)}(t_k) \cdot \mathbf{l}^{(i)}(t_k)}{d^{(i)}(t_k)} + \delta_{t,i}$$

其中 $\delta_{t,i}$ 是第 $i$ 颗卫星的等效钟漂（单位 Hz，视为常数）。

### 9.2 未知数扩展

$$\mathbf{x} = [x_u, y_u, z_u, \delta_{t,1}, \delta_{t,2}, \ldots, \delta_{t,N}]^T$$

- 位置：3 个
- 钟漂：$N$ 个（每颗卫星 1 个）
- 总计：$3 + N$ 个

### 9.3 扩展 Jacobian

$$\mathbf{H}_{\text{ext}} = \begin{bmatrix}
\frac{\partial f_D^{(1)}(t_1)}{\partial x_u} & \frac{\partial f_D^{(1)}(t_1)}{\partial y_u} & \frac{\partial f_D^{(1)}(t_1)}{\partial z_u} & 1 & 0 & \cdots & 0 \\[4pt]
\frac{\partial f_D^{(1)}(t_2)}{\partial x_u} & \frac{\partial f_D^{(1)}(t_2)}{\partial y_u} & \frac{\partial f_D^{(1)}(t_2)}{\partial z_u} & 1 & 0 & \cdots & 0 \\[4pt]
\vdots & \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\[4pt]
\frac{\partial f_D^{(N)}(t_K)}{\partial x_u} & \frac{\partial f_D^{(N)}(t_K)}{\partial y_u} & \frac{\partial f_D^{(N)}(t_K)}{\partial z_u} & 0 & 0 & \cdots & 1
\end{bmatrix}_{M \times (3+N)}$$

因为 $\partial f_D/\partial \delta_{t,i} = 1$（钟漂直接加到多普勒上），所以钟漂对应的 Jacobian 列是**稀疏的 0/1 矩阵**：第 $i$ 颗卫星的观测行在第 $(3+i)$ 列为 1，其余为 0。

---

## 10. 与论文的对应

| 论文 (Shahcheraghi & Kassas, 2024) | 本文 |
|------------------------------------|------|
| Acquisition: Sequential MSD + FFT-based GLR | $\leftarrow$ 提供 $f_D$ 测量值 |
| Doppler Tracking: 3rd-order KF | $\leftarrow$ 提纯 $f_D$ 测量值 |
| **Positioning: NLS filter** | **$\leftarrow$ 本文核心（Step 3）** |
| State vector: position + clock drifts | $\mathbf{x} = [x_u,y_u,z_u,\delta_{t,1}\ldots\delta_{t,N}]$ |
| Initial error: 179 km | 本文: 50 km |
| Iterations: 7 | 取决于初始偏差和几何 |
| Final error: 9.52 m (2D) | 取决于测量噪声和卫星数 |
| Stopping threshold: $10^{-9}$ | $\|\Delta\mathbf{p}\| < 10^{-6}$ m |

---

## 参考文献

1. Shahcheraghi, S. & Kassas, Z.M. (2024). "A Computationally Efficient Approach for Acquisition and Doppler Tracking for PNT With LEO Megaconstellations." *IEEE Signal Processing Letters*, 31, 2400–2404.
2. Neinavaie, M., Khalife, J. & Kassas, Z.M. (2022). "Acquisition, Doppler Tracking, and Positioning with Starlink LEO Satellites: First Results." *IEEE Trans. Aerosp. Electron. Syst.*, 58(3), 2606–2610.
