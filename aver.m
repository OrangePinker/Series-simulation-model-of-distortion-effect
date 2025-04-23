% 读取Excel文件
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\distortion.xlsx'; % 替换为您的文件名
data = readmatrix(filename);

% 提取电流值和光谱数据
currents = data(1, :);          % 第一行为电流值
spectra = data(2:end, :);       % 第2-513行为光谱数据

% 获取唯一电流值（保持原始顺序）
[unique_currents, ~, ic] = unique(currents, 'stable');

% 预分配存储空间
avg_spectra = zeros(size(spectra, 1), length(unique_currents));

% 计算每个电流值的平均光谱
for k = 1:length(unique_currents)
    cols = (currents == unique_currents(k));       % 查找当前电流对应的列
    avg_spectra(:, k) = mean(spectra(:, cols), 2); % 计算列平均
end

% 对平均光谱进行四舍五入取整
avg_spectra = round(avg_spectra);  % 新增的取整操作

% 构建输出矩阵
output = [unique_currents; avg_spectra];

% 写入结果到新Excel文件
writematrix(output, 'averaged_result.xlsx');

%%
% 读取Excel数据，范围根据实际情况调整
data = readmatrix('normal.xlsx', 'Range', 'A2:T513'); % 20列，512行

% 计算每行的极差（最大值 - 最小值）
R = max(data, [], 2) - min(data, [], 2);

% 设置极差法系数d（n=20时，d≈3.735）
d = 3.735;

% 计算标准差估计
sigma = R / d;

% 显示结果
disp('每个测量点的误差估计：');
disp(sigma);

% 可选：保存结果到Excel
writematrix(sigma, 'error_estimates.xlsx');

%%


