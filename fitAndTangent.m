function [fitParams, tangentParams] = fitAndTangent(filename, x0)
    % 读取Excel文件
    data = xlsread(filename);
    
    % 获取x坐标和y坐标
    x_data = data(1, :);    % 第1行作为x坐标
    y_data = data(514, :);  % 第514行作为y坐标
    
    % 进行三次多项式拟合
    p = polyfit(x_data, y_data, 3);
    
    % 计算拟合值
    y_fit = polyval(p, x_data);
    
    % 计算切点y0
    y0 = polyval(p, x0);
    
    % 计算切线斜率m
    m = 3*p(1)*x0^2 + 2*p(2)*x0 + p(3);
    
    % 生成切线的x值
    x_tangent = linspace(x0, x0 + 400, 400);
    
    % 计算切线的y值
    y_tangent = m*(x_tangent - x0) + y0;
    
%     % 绘图
%     figure;
%     plot(x_data, y_data, 'r--', 'LineWidth', 2); % 原始数据
%     hold on;
%     plot(x_data, y_fit, 'm-', 'LineWidth', 1); % 拟合曲线
%     plot(x_tangent, y_tangent, 'b-.', 'LineWidth', 1.5); % 切线
%     plot(x0, y0, 'dk', 'MarkerSize', 8); % 切点
%     grid on;
%     xlabel('X-ray Tube Current [μA]');
%     ylabel('Photon Counts');
%     legend('Original Data', 'Cubic Fit', 'Tangent Line', 'Tangent Point', 'Location', 'Best');
%     xlim([150 650]); % 设置x轴范围
%     ylim([1410000 2000000]); % 设置y轴范围

    % 输出拟合方程和切线方程的参数
    fprintf('拟合方程: y = %.4fx^3 + %.4fx^2 + %.4fx + %.4f\n', p(1), p(2), p(3), p(4));
    fprintf('切线方程: y = %.4f x + %.4f\n', m, y0 - m*x0);

    % 返回拟合参数和切线参数
    fitParams = p;
    tangentParams = [m, y0 - m*x0];
end