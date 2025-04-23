function [params1, params2] = fitTwoWeibull(filename)
% 读取Excel文件
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

% % 绘制数据的概率密度折线图
% figure;
% plot(values, prob_density, 'r-', 'DisplayName', 'Original Data','LineWidth', 1.5);
% xlabel('Value [1-512]');
% ylabel('Probability Density');
% hold on;

% 选择要检查的分布类型'Normal', 
distributions = {'Exponential', 'Gamma', 'Lognormal', 'Weibull'}; % 可以添加更多分布

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

%         % 绘制拟合结果
%         plot(x_values, y_values, 'DisplayName', dist);

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

% hold off;
% legend show;

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
% 输出参数
params1 = paramValues(1);
params2 = paramValues(2);
end