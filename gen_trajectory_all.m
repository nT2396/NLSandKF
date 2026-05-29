%% generate_tle_mat.m — 从TLE生成timetable格式的.mat星历文件 (全部卫星)
% 输出: timetableSatelliteTrajectory.mat (positionTT, velocityTT)

clear; clc;

%% ===== 参数配置 =====
startTime  = datetime(2026, 5, 25, 5, 50, 0);
stopTime   = datetime(2026, 5, 25, 6, 0, 0);
sampleTime = 10;
tleFile    = 'Sat_tle.txt';
outFile    = 'timetableSatelliteTrajectory_300.mat';

%% ===== Step 1: 创建卫星场景, 加载 TLE =====
fprintf('=== 加载 TLE 文件: %s ===\n', tleFile);
sc   = satelliteScenario(startTime, stopTime, sampleTime);
sats = satellite(sc, tleFile);
N_sat = numel(sats);
fprintf('成功加载 %d 颗卫星\n', N_sat);

%% ===== Step 2: 计算各卫星 ECEF 位置和速度 =====
t   = startTime : seconds(sampleTime) : stopTime;
N_t = numel(t);
fprintf('\n时间范围: %s → %s\n', string(t(1)), string(t(end)));
fprintf('采样间隔: %d 秒, 共 %d 个时间点\n', sampleTime, N_t);

posCell = cell(N_t, N_sat);
velCell = cell(N_t, N_sat);

fprintf('\n计算卫星轨道 (ECEF)...\n');
for i = 1:N_sat
    p_mat = zeros(N_t, 3);
    v_mat = zeros(N_t, 3);
    for k = 1:N_t
        [p, v] = states(sats(i), t(k), 'CoordinateFrame', 'ecef');
        p_mat(k, :) = p';
        v_mat(k, :) = v';
    end
    posCell(:, i) = num2cell(p_mat, 2);
    velCell(:, i) = num2cell(v_mat, 2);
    if mod(i, 10) == 0 || i == N_sat
        fprintf('  卫星 %d / %d\n', i, N_sat);
    end
end

%% ===== Step 3: 构建 timetable =====
timeDuration = seconds(0 : sampleTime : (N_t-1)*sampleTime)';
vNames = compose('Sat%d', 1:N_sat);

positionTT = table2timetable(cell2table(posCell, 'VariableNames', vNames), ...
                             'RowTimes', timeDuration);
velocityTT = table2timetable(cell2table(velCell, 'VariableNames', vNames), ...
                             'RowTimes', timeDuration);

%% ===== Step 4: 保存 =====
save(outFile, 'positionTT', 'velocityTT');
fprintf('\n已保存: %s (%d 颗卫星, %d 时间点)\n', outFile, N_sat, N_t);
load("timetableSatelliteTrajectory_4Sat.mat","positionTT","velocityTT");
sat = satellite(sc,positionTT,velocityTT,"CoordinateFrame","ecef");
play(sc);