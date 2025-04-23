% 读取Excel文件（请替换为您的文件路径）
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\fsc250331\fsc250331\paper\PAPER.xlsx';       % Excel文件名

sheet = 'Sheet1';                  % 工作表名称
rowNumber = 514;                   % 要读取的行号

% 读取指定行数据（建议先检查数据是否存在标题行）
fullData = readmatrix(filename, 'Sheet', sheet);
selectedRow = fullData(rowNumber, :);

% 定义全新x轴坐标 (500到600步长20)
x = 500:20:600;  % 自动生成 [500,520,540,560,580,600]

% 数据列分组验证
numGroups = floor(length(selectedRow)/6);
assert(numGroups >=1, '数据列数不足以形成至少1个完整分组');
colors = lines(numGroups);  % 使用内置色谱

% 创建图形窗口
figure;
hold on;
% 预定义图例名称（新增部分）
legendNames = {'10 sheets', '40 sheets', '70 sheets', '100 sheets'};
% 每组数据绑定新x坐标
for i = 1:numGroups
    % 定位列范围
    cols = (i-1)*6 + (1:6);
    
    % 提取对应数据段
    plotData = selectedRow(cols);
    
    % 绘制带标签的曲线
    plot(x, plotData,...
        'Color', colors(i,:),...
        'LineWidth', 1.5,...
        'Marker', 'o',...
        'MarkerFaceColor', 'white',...
        'DisplayName', sprintf('Thickness %d', i)); 
end

% 图形美化设置
hold off;
xlabel('Parameter Value (Unit)', 'FontWeight','bold'); % 请替换实际单位和意义
ylabel('Measurement Value (Unit)', 'FontWeight','bold'); 
title('Thickness experiment', 'FontSize',12);
legend('Location', 'bestoutside', 'NumColumns',ceil(numGroups/4));
grid on;

% 精确设置x轴显示
set(gca,...
    'XLim', [480 620],...          % 留出边距的显示范围
    'XTick', 500:20:600,...        % 精确刻度步长
    'XMinorTick', 'on',...         % 显示次级刻度
    'YLim', [1720000, 1880000],... % 自动y轴适应
    'FontName', 'Arial');  

% 给折线图添加x轴和y轴标签
xlabel('Tube Current');  % 在这里给 x 轴添加标签
ylabel('Photon Counts');

% 预定义图例名称（新增部分）
legendNames = {'10 sheets', '40 sheets', '70 sheets', '100 sheets'};

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

%%
% 导出图形为 JPG 图片
print('-djpeg', '纸张厚度实验', '-r0');

% 关闭图形
close(figureHandle);