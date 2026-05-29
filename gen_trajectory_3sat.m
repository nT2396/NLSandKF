%% gen_tle_v3.m — 从STARLINK_2中选取3颗卫星, 生成轨迹
% 选取: Sat59(STARLINK-1270), Sat74(STARLINK-1274), Sat86(STARLINK-1309)
% 时间: 2026-05-25 08:38:00 ~ 08:42:00 UTC
% 输出: traj_3Sat.mat

clear all

%% ===== 参数配置 =====
startTime  = datetime(2026, 5, 25, 8, 38, 0, 'TimeZone', 'UTC');
stopTime   = datetime(2026, 5, 25, 8, 42, 0, 'TimeZone', 'UTC');
sampleTime = 0.2;
tleFile    = 'STARLINK_2.txt';
outFile    = 'traj_3Sat_step02.mat';
selSat     = [59, 74, 86];  % STARLINK-1270, STARLINK-1274, STARLINK-1309

%% ===== Step 1: 加载全部 TLE, 筛选指定卫星 =====
fprintf('=== 加载 TLE 文件: %s ===\n', tleFile);
sc   = satelliteScenario(startTime, stopTime, sampleTime);
sats = satellite(sc, tleFile);
N_all = numel(sats);
fprintf('TLE中共 %d 颗卫星, 选取第 %s 号\n', N_all, num2str(selSat));

sats = sats(selSat);
N_sat = numel(sats);
fprintf('选取 %d 颗卫星:\n', N_sat);
for i = 1:N_sat
    fprintf('  [Sat%d] %s\n', selSat(i), sats(i).Name);
end

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
    fprintf('  Sat%d / %d : %s (%d 个点)\n', selSat(i), N_sat, sats(i).Name, N_t);
end

%% ===== Step 3: 构建 timetable =====
timeDuration = seconds(0 : sampleTime : (N_t-1)*sampleTime)';
vNames = compose('Sat%d', selSat);

positionTT = table2timetable(cell2table(posCell, 'VariableNames', vNames), ...
                             'RowTimes', timeDuration);
velocityTT = table2timetable(cell2table(velCell, 'VariableNames', vNames), ...
                             'RowTimes', timeDuration);

%% ===== Step 4: 验证 =====
fprintf('\n=== 格式验证 ===\n');
fprintf('positionTT: %d x %d, 变量: ', height(positionTT), width(positionTT));
for i = 1:N_sat, fprintf('%s ', vNames{i}); end
fprintf('\n');
val = positionTT{1, 1};
fprintf('数据格式: class=%s, size=[%s]\n', class(val), num2str(size(val)));

%% ===== Step 5: 保存 =====
save(outFile, 'positionTT', 'velocityTT');
fprintf('\n已保存: %s (%d 颗卫星, %d 时间点)\n', outFile, N_sat, N_t);

%% ===== Step 6: 3D 可视化 =====
% sc = satelliteScenario;
% load(outFile, 'positionTT', 'velocityTT');
% sat = satellite(sc, positionTT, velocityTT, 'CoordinateFrame', 'ecef');
% play(sc);
