%%
%拟合函数的测试
Total_filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\total.xlsx';
x0 = 200;
[fitParams, tangentParams] = fitAndTangent(Total_filename, x0);
% 提取参数
fitParams_a = fitParams(1);
fitParams_b = fitParams(2);
fitParams_c = fitParams(3);
fitParams_d = fitParams(4);
tangentParams_m = tangentParams(1);
tangentParams_e = tangentParams(2);
%%
[params1, params2] = fitTwoWeibull('C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx');
a = params1
b = params2

%%
clc; clear;
% 导入Excel文件
%拟合光子数曲线
filename = 'total.xlsx';
data = xlsread(filename);

% 获取x坐标和y坐标
y_data = data(514, :);  % 第514行作为y坐标
x_data = data(1, :);    % 第1行作为x坐标

% 确保数据是列向量并且没有NaN或Inf
x_data = x_data(:);
y_data = y_data(:);
validIdx = isfinite(x_data) & isfinite(y_data);
x_data = x_data(validIdx);
y_data = y_data(validIdx);

% 数据标准化
x_data = (x_data - mean(x_data)) / std(x_data);

% 定义要尝试的degree范围
degrees = 1:5;

% 初始化R_squared数组和存储拟合系数
R_squared = zeros(size(degrees));
coefficients = cell(1, length(degrees));

% 循环进行拟合和计算R²
for i = 1:length(degrees)
    degree = degrees(i);
    p = polyfit(x_data, y_data, degree);
    coefficients{i} = p;
    y_fit = polyval(p, x_data);
    SSR = sum((y_fit - mean(y_data)).^2);
    SST = sum((y_data - mean(y_data)).^2);
    R_squared(i) = SSR / SST;
    fprintf('Degree %d: R² = %f\n', degree, R_squared(i));
end

% 输出拟合方程
for i = 1:length(degrees)
    degree = degrees(i);
    p = coefficients{i};
    fprintf('Degree %d 拟合方程：y = ', degree);
    for j = 1:length(p)
        fprintf('%f', p(j));
        if j < length(p)
            fprintf('x^%d + ', length(p)-j);
        end
    end
    fprintf('\n');
end

% 可视化所有拟合结果
figure;
plot(x_data, y_data, 'o', 'DisplayName', '原始数据');
hold on;
colors = ['r','g','b','c','m'];
for i = 1:length(degrees)
    degree = degrees(i);
    p = coefficients{i};
    y_fit = polyval(p, x_data);
    plot(x_data, y_fit, '-', 'Color', colors(i), 'DisplayName', ['Degree ', num2str(degree)]);
end
legend show;
xlabel('x');
ylabel('y');
title('不同阶数多项式拟合的比较');

%%
%拟合光子数和电流增长一次直线
% 导入Excel文件
clc; clear;

filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\total.xlsx';
data = xlsread(filename);

% 获取数据的行数和列数
[numRows, numCols] = size(data);

% 获取x坐标和y坐标
y_data = data(514, :);  % 第514行作为y坐标
x_data = data(1, :);    % 第1行作为x坐标

% 绘制折线图
figure; % 新增图形
plot(x_data, y_data, 'm-', 'LineWidth', 1); % 将线条加粗
grid on;
hold on;

% 进行三次多项式拟合
p = polyfit(x_data, y_data, 3); 
y_fit = polyval(p, x_data); % 计算拟合值

% 选择切点 x0
x0 = 200;
y0 = polyval(p, x0); % 计算切点 y0
m = 3*p(1)*x0^2 + 2*p(2)*x0 + p(3); % 计算切线斜率
x_tangent = linspace(x0, x0 + 60, 100); % 生成切线的 x 值
y_tangent = m*(x_tangent - x0) + y0;

% 绘制拟合线
plot(x_data, y_fit, 'r--', 'LineWidth', 2); % 拟合线
plot(x_tangent, y_tangent, 'b-.', 'LineWidth', 1.5); % 切线
plot(x0, y0, 'dk', 'MarkerSize', 8); % 标出切点
legend('Tangent Equation', 'Quadratic Fitting', 'Original Data', 'Location', 'Best'); % 添加图例

% 给折线图添加x轴和y轴标签
xlabel('X-ray Tube Current [μA]');  % 在这里给 x 轴添加标签
ylabel('Photon Counts');

% 输出拟合方程
a = p(1);
b = p(2);
c = p(3);
d = p(4);
fprintf('拟合方程: y = %.4fx^3 + %.4fx^2 + %.4fx + %.4f\n', a, b, c, d);
% 计算新的截距 new_c
new_c = -m*x0 + y0;
fprintf('切线方程：y = %.4f x + %.4f\n', m, new_c);

% 设置坐标轴属性
set(gca, 'Box', 'off', ...
         'Layer', 'top', ...
         'LineWidth', 1, ...
         'XGrid', 'off', 'YGrid', 'off', ...
         'TickDir', 'out', 'TickLength', [0.01 0.01], ...
         'XMinorTick', 'off', 'YMinorTick', 'off', ...
         'XColor', [.1 .1 .1], 'YColor', [.1 .1 .1]);

% 获取 x,y 轴标签的句柄
hXLabel = get(gca, 'XLabel');
hYLabel = get(gca, 'YLabel');
hTitle = get(gca, 'Title');

% 设置字体和字号
set(gca, 'FontName', 'Arial', 'FontSize', 10);
set(hTitle, 'FontSize', 12, 'FontWeight', 'bold');

% 设置x轴和y轴的显示范围
xlim([150 300]); % 设置x轴范围
ylim([1418000 1750000]); % 设置y轴范围

% 设置背景颜色
set(gcf, 'Color', [1 1 1]);

% 添加上、右框线
xc = get(gca, 'XColor');
yc = get(gca, 'YColor');
unit = get(gca, 'units');
ax = axes('Units', unit, ...
           'Position', get(gca, 'Position'), ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...
           'XColor', xc, ...
           'YColor', yc);
set(ax, 'linewidth', 1, ...
        'XTick', [], ...
        'YTick', []);

% 获取当前图形的句柄
figureHandle = gcf;

% 设置图形的宽度和高度
set(figureHandle, 'PaperUnits', 'points');
set(figureHandle, 'PaperPosition', [0 0 600 400]);
% 导出图形为 JPG 图片
print('-djpeg', '理论电流-光子数关系', '-r0');

% 关闭图形
close(figureHandle);
%%
%等比放大
clc; clear;

% 导入Excel文件
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx';
data = xlsread(filename);

% 获取数据的行数，假设只有一列数据
[numRows, numCols] = size(data);

% 确保数据行数足够
if numRows < 513
    error('数据行数不足，至少应包含 513 行。');
end

% 获取原始数据：第2行到第513行
original_data = data(2:513, 1);  % 取第一列的第2到513行数据

% 将原始数据分到512个数组内
num_bins = 512;
bins = cell(1, num_bins);

% 将数据分到512个bin中
for i = 1:num_bins
    if i <= length(original_data)
        bins{i} = original_data(i);  % 当前bin取当前值
    else
        bins{i} = 0;  % 对于不足的bin，填充为 0
    end
end

% 计算增量
initial_current = 200;
Current_Value = 280;
incre = (5013.5214 * Current_Value + 421054.5697) / (5013.5214 * initial_current + 421054.5697);

% 将导入的数据依次与incre相乘
scaled_data = zeros(1, num_bins);  % 创建存储结果的数组
for i = 1:num_bins
    scaled_data(i) = bins{i} * incre;  % 乘以incre
end

% 绘制结果
figure; % 创建一个新图形窗口
hold on; % 保持当前图形，以便绘制多个数据系列
plot(original_data(1:num_bins), 'r-', 'LineWidth', 1.5); % 绘制原始数据，红色线
plot(scaled_data, 'b-', 'LineWidth', 1.5); % 绘制缩放后的数据，蓝色线
hold off; % 释放图形保持

% 添加图例
legend('Original Data', 'Scaled Data');

% 添加标签和标题
xlabel('Bin Index'); % x轴标签
ylabel('Data Value'); % y轴标签
title('Original Data and Scaled Data vs. Bin Index'); % 图形标题
grid on; % 添加网格

%%
%分布拟合
clc; clear;

% 导入Excel文件
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx';
data = xlsread(filename);

% 获取数据的行数
[numRows, numCols] = size(data);

% 确保数据行数足够
if numRows < 513
    error('数据行数不足，至少应包含 513 行。');
end

% 获取 Y 数组：第 2 行到第 513 行，第一列
Y = data(2:513, 1);

% 假设 Y 的索引代表值，元素代表频数
values = 0:length(Y)-1; % 值从0到511
counts = Y; % 频数

% 生成数据向量
data_vector = repelem(values, counts);

% 检查并过滤掉非正值的数据
positive_data = data_vector(data_vector > 0);

% 确保positive_data是数值型的列向量，并且不含NaN或Inf
if ~isnumeric(positive_data) || any(isnan(positive_data)) || any(isinf(positive_data))
    error('positive_data包含非数值或无效值。');
end
positive_data = positive_data(:);

% 通知用户零值的数量
zero_counts = sum(counts == 0);
fprintf('数据中存在 %d 个零值\n', zero_counts);

% 必要时处理零值，可以选择跳过或替换
if isempty(positive_data)
    error('数据中不包含正值，无法进行分布拟合。');
end

% 计算概率密度
total = sum(counts);
prob_density = counts / total;

% 绘制数据的概率密度折线图
figure;
plot(values, prob_density, 'r-', 'DisplayName', 'Original Data','LineWidth', 1.5);
xlabel('Value [1-512]');
ylabel('Probability Density');
hold on;

% 选择要检查的分布类型
distributions = {'Normal', 'Exponential', 'Gamma', 'Lognormal', 'Weibull'}; % 可以添加更多分布

% 存储拟合优度信息
fitResults = struct();
aicValues = zeros(length(distributions), 1); % 存储每种分布的AIC值
fitSuccess = true(length(distributions), 1); % 记录每种分布拟合是否成功

% 对每种分布进行拟合并计算AIC
for i = 1:length(distributions)
    dist = distributions{i};
    
    try
        pd = fitdist(positive_data, dist);
        
        % 生成拟合的概率密度函数
        x_values = linspace(min(positive_data), max(positive_data), 100);
        y_values = pdf(pd, x_values);
        
        % 绘制拟合结果
        plot(x_values, y_values, 'DisplayName', dist);
        
        % 计算AIC（赤池信息量）
        logLikelihood = sum(log(pdf(pd, positive_data)));
        numParams = numel(pd.ParameterNames);
        AIC = -2 * logLikelihood + 2 * numParams;
        
        % 存储结果
        fitResults.(dist) = struct('pdf', y_values, 'AIC', AIC, 'params', pd);
        aicValues(i) = AIC; % 存储AIC值
    catch ME
        fprintf('拟合 %s 时出错: %s\n', dist, ME.message);
        fitSuccess(i) = false;  % 记录拟合失败
        aicValues(i) = NaN;  % 无效的AIC值
    end
end

hold off;
legend show;

% 比较拟合优度，找到AIC最小的分布
validAIC = aicValues(~isnan(aicValues));
if ~isempty(validAIC)
    [~, minIndex] = nanmin(aicValues);
    bestDistribution = distributions{minIndex};
else
    bestDistribution = '无合适的分布';
end

% 输出结果，包括各个分布的 AIC 值
fprintf('各个分布的 AIC 值:\n');
for i = 1:length(distributions)
    if fitSuccess(i)
        fprintf('%s: AIC = %.4f\n', distributions{i}, aicValues(i));
    else
        fprintf('%s: 拟合失败\n', distributions{i});
    end
end

% 输出最接近的分布及其参数
if ~isempty(validAIC)
    pd_best = fitResults.(bestDistribution).params;
    paramNames = pd_best.ParameterNames;
    paramValues = pd_best.ParameterValues;
    fprintf('最佳分布的参数:\n');
    for i = 1:length(paramNames)
        fprintf('%s: %.4f\n', paramNames{i}, paramValues(i));
    end
    fprintf('数据最接近的分布是: %s（AIC = %.4f）\n', bestDistribution, aicValues(minIndex));
else
    fprintf('所有分布拟合失败，无法确定最优分布。\n');
end


%%
% 脉冲堆积
clc; clear;

% 记录运行开始时间
tic;

% 导入Excel文件
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx';
data = xlsread(filename);

% 获取数据的行数，假设只有一列数据
[numRows, numCols] = size(data);

% 确保数据行数足够
if numRows < 513
    error('数据行数不足，至少应包含 513 行。');
end

% 获取原始数据：第2行到第513行
original_data = data(2:513, 1);  % 取第一列的第2到513行数据

% 将原始数据分到512个数组内
num_bins = 512;
bins = cell(1, num_bins);

% 将数据分到512个bin中
for i = 1:num_bins
    if i <= length(original_data)
        bins{i} = original_data(i);  % 当前bin取当前值
    else
        bins{i} = 0;  % 对于不足的bin，填充为 0
    end
end

% 计算增量
initial_current = 200;
Current_Value = 600;
incre = (5013.5214 * Current_Value + 421054.5697) / (5013.5214 * initial_current + 421054.5697);

% 将导入的数据依次与incre相乘
scaled_data = zeros(1, num_bins);  % 创建存储结果的数组
for i = 1:num_bins
    scaled_data(i) = bins{i} * incre;  % 乘以incre
end

% 计算原始总数
scaled_data = round(scaled_data);
original_total = sum(scaled_data);
incre_data = zeros(1, num_bins);
for i = 1:num_bins
    incre_data(i) = bins{i} * incre;  % 乘以incre
end

% 初始化进度条
processed_counts = 0;
percent_step = 1;
last_update_percent = 0;
h = waitbar(0,'Processing...');

% 新数组初始化
new_array = zeros(1, 512);

% 韦布尔分布参数
scale_param = 181.3227;
shape_param = 2.3587;

% 定义每次取出数量的概率分布
% 例如，90%的概率取1，10%的概率取2
Probability = (0.0146 * Current_Value^3 + -25.3121 * Current_Value^2 + 13316.2956 * Current_Value^1 - 342061.0848)/(5013.5214 * Current_Value + 421054.5697)
prob_take_1 = Probability;
prob_take_2 = Probability+(1 - Probability)*0.81;
prob_take_3 = 1.00;

% 循环直到scaled_data中的所有数被取完
while any(scaled_data > 0)

    % 更新processed_counts
    zeroCount = sum(scaled_data == 0);

    if toc>60
        break;
    end

    % 更新进度条
    current_percent = (zeroCount / 512) * 100;
    if current_percent >= last_update_percent + percent_step
        waitbar(current_percent / 100, h, sprintf('Processing: %.2f%%', current_percent));
        last_update_percent = current_percent;
    end

    % 生成韦布尔分布的随机数来决定取出哪个数
    % 需要将随机数映射到1到512之间
    random_num = wblrnd(scale_param, shape_param, 1, 1);
    bin_index = round(random_num);

    if bin_index < 7 || bin_index > 512 || scaled_data(bin_index) < 1
        continue;
    end


    % 决定取出的数量
    r = rand;
    if r < prob_take_1
        take_num = 1;
    elseif r >= prob_take_1 && r < prob_take_2
        take_num = 1.5; 

    else
        take_num = 2;

    end
    
    while (take_num > 1)

        random_num2 = wblrnd(scale_param, shape_param, 1, 1);
        bin_index2 = round(random_num2);

        if bin_index2 < 1 || bin_index2 > 512 || scaled_data(bin_index2) < 1
            continue;  % 跳过此次循环，不进行赋值
        end

        break;

    end

    % 检查取出的数量是否大于1
    if take_num == 1
        % 取出一个数
        if scaled_data(bin_index) > 0
            new_array(bin_index) = new_array(bin_index) + 1;
            scaled_data(bin_index) = scaled_data(bin_index) - 1;
        end

    elseif take_num == 1.5
        compare_bin = max(bin_index2,bin_index);
        if compare_bin < 7
            continue;  % 跳过此次循环，不进行赋值
        end
        if scaled_data(compare_bin) > 0
            new_array(bin_index) = new_array(bin_index) + 1;
            scaled_data(bin_index) = scaled_data(bin_index) - 1;
            scaled_data(bin_index2) = scaled_data(bin_index2) - 1;
        end
        
    else
        % 随机选择take_num个数
        selected_indices = bin_index2 + bin_index;
        if selected_indices < 7
            continue;  % 跳过此次循环，不进行赋值
        end
        %sum_indices = sum(selected_indices);
        if selected_indices > 512
            scaled_data(bin_index2) = scaled_data(bin_index2) - 1;
            scaled_data(bin_index) = scaled_data(bin_index) - 1;
            continue;  % 跳过此次循环，不进行赋值
        end
        if selected_indices <= 512
            % 检查这些bin是否有足够的数量
            sufficient = true;
            for i = 1:take_num
                if scaled_data(bin_index2) < 1 || scaled_data(bin_index) < 1
                    sufficient = false;
                    break;
                end
            end
            if sufficient
                new_array(selected_indices) = new_array(selected_indices) + 1;
                scaled_data(bin_index2) = scaled_data(bin_index2) - 1;
                scaled_data(bin_index) = scaled_data(bin_index) - 1;
                
            end
        end
    end
end

% 计算新总数
new_total = sum(new_array);

% 关闭进度条
close(h);

% 输出总数
fprintf('Original Total: %f\n', original_total);
fprintf('New Total: %f\n', new_total);

% 输出运行时间
fprintf('Total runtime: %f seconds.\n', toc);

% 绘制结果
figure; % 创建一个新图形窗口
hold on; % 保持当前图形，以便绘制多个数据系列
plot(incre_data, 'r-', 'LineWidth', 2); % 绘制原始数据，红色线
plot(new_array, 'b-', 'LineWidth', 2); % 绘制新数组，蓝色线
hold off; % 释放图形保持

% 添加图例
legend('Original Data', 'Pile up Data');

% 添加标签和标题
xlabel('Channel [1-512]'); % x轴标签
ylabel('Photon Counts'); % y轴标签

grid on; % 添加网格

% 设置整体属性
set(gcf, 'Color', [1 1 1]);
% 获取当前图形的句柄
figureHandle = gcf;

% 设置图形的宽度和高度
set(figureHandle, 'PaperUnits', 'points');
set(figureHandle, 'PaperPosition', [0 0 600 400]);
% % 导出图形为 JPG 图片
% print('-djpeg', '脉冲堆积0113', '-r0');
% 
% % 关闭图形
% close(figureHandle);

