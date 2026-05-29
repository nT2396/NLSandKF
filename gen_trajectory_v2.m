%% gen_tle_mat_v2.m — 从93颗中选取指定4颗卫星, 生成 .mat
% 选取: Sat70(STARLINK-1319), Sat71(STARLINK-1266), Sat82(STARLINK-1263), Sat88(STARLINK-1286)
% 输出: timetableSatelliteTrajectory_4Sat.mat

clear all

%% ===== 参数配置 =====
startTime  = datetime(2026, 5, 25, 10, 00, 0,'TimeZone','UTC');
stopTime   = datetime(2026, 5, 25, 10, 10, 0,'TimeZone','UTC');
sampleTime = 10;
tleFile    = 'STARLINK_TLE_part4.txt';
outFile    = 'timetableSatelliteTrajectory_part4.mat';
selSat     = [70, 71, 82, 88];  % 选取这4颗卫星

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

% sc = satelliteScenario;
load("timetableSatelliteTrajectory_part4.mat","positionTT","velocityTT");
sat = satellite(sc,positionTT,velocityTT,"CoordinateFrame","ecef");
play(sc);