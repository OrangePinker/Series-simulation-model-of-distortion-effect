% 极化效应

% 导入Excel文件
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx';
data = xlsread(filename);

% 记录运行开始时间
tic;

% 获取数据的行数，假设只有一列数据
[numRows, numCols] = size(data);

% 确保数据行数足够
if numRows < 513
    error('数据行数不足，至少应包含 513 行。');
end

% 获取原始数据：第2行到第513行
original_data = data(2:513, 1);  % 取第一列的第2到513行数据

% 参数设置
Current_Value = 200;
Upper_Limit = 600;
Lower_Limit = 200;
Polarization_Q0 = (8.84e-6 - 8.6e-6)* (Current_Value - Lower_Limit) / (Upper_Limit - Lower_Limit)+8.6e-6
t = calculate_time(Polarization_Q0, 1e-6, 500, 0.01);
K = 3e6;
Attenuation_Rate = 1 - K * t

% 创建新数组
pola_data = zeros(1, 512);

% 初始化变量
pending_data = 0;
last_non_zero_index = 0;

% 处理每个原始数据点
for i = 1:512
    new_index = i * Attenuation_Rate;
    
    % 舍弃新索引小于7的数据
    if new_index < 7
        continue;  % 跳过此次循环，不进行赋值
    end
    
    new_index = round(new_index);
    
    if pola_data(new_index) ~= 0
        if pending_data ~= 0
            % 分配pending_data到last_non_zero_index和new_index之间的位置
            start = last_non_zero_index + 1;
            end_pos = new_index;
            if start <= end_pos
                num_positions = end_pos - start + 1;
                allocation = pending_data / num_positions;
                for k = start:end_pos
                    pola_data(k) = pola_data(k) + allocation;
                end
                % 调试信息
                %disp(['分配 pending_data: ', num2str(pending_data), ' 到位置 ', num2str(start), ' 到 ', num2str(end_pos)]);
            end
            pending_data = 0;  % 分配完后清零
        end
        % 保留new_data(new_index)，跳过放入original_data(i)
        pending_data = original_data(i);  % 存入新的original_data(i)
        last_non_zero_index = new_index;
    else
        pola_data(new_index) = pola_data(new_index) + original_data(i);

    end
end

% 绘制图形
figure;
plot(1:512, original_data, 'r-', 'DisplayName', 'Original Data','LineWidth', 2);
hold on;
plot(1:512, pola_data, 'b-', 'DisplayName', 'Polarization Data','LineWidth', 2);
hold off;
legend show;
xlabel('Channel [1-512]');
ylabel('Photon Counts');

grid on;

% 计算并输出数组的和
sum_original = sum(original_data);
sum_new = sum(pola_data);

disp(['原始数据的和: ', num2str(sum_original)]);
disp(['新数据的和: ', num2str(sum_new)]);

% 输出运行时间
fprintf('Total runtime: %f seconds.\n', toc);
% % 将新数据导出到Excel的第2到第513行，第一列
% writematrix(new_data.', 'output.xlsx', 'Range', 'A2');

%%
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
% 定义 Polarization_Q0 的范围
Polarization_Q0 = 8e-6 : 1e-11 : 8.85399e-6;

% 初始化 t 向量
t = zeros(size(Polarization_Q0));

% 循环计算每个 Polarization_Q0 对应的 t
for i = 1:length(Polarization_Q0)
    t(i) = calculate_time(Polarization_Q0(i), 1e-6, 500, 0.01);
end

% 绘制 Polarization_Q0 和 t 的关系图
figure;
plot(Polarization_Q0, t, 'r-','LineWidth', 1.5);
xlabel('Polarization Q0');
ylabel('Time Of Flight');
grid on;
