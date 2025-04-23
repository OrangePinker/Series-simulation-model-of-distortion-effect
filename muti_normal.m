Total_filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\fsc250331\fsc250331\PLA\PLA.xlsx';
x0 = 200;
[fitParams, tangentParams] = fitAndTangent(Total_filename, x0);
tangentParams_m = tangentParams(1);
tangentParams_e = tangentParams(2);

% filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\normal.xlsx';
filename = Total_filename ;

% 导入并验证数据
data = xlsread(filename);
[numRows, numCols] = size(data);


% 参数设置
initial_current = 200;
current_values = 200:10:600;  % 生成40个电流值（200到600）
total_columns = 1 * length(current_values);  % 总列数800
output_matrix = zeros(513, total_columns);    % 预分配内存

% 双层循环处理
output_col = 1;  % 输出列计数器
for orig_col = 1:1          % 原始数据列循环
    raw_data = data(2:513, orig_col);  % 获取原始光谱数据
    
    for cv_idx = 1:length(current_values)  % 电流值循环
        Current_Value = current_values(cv_idx);
        
        % 计算增量系数
        incre = (tangentParams_m * Current_Value + tangentParams_e) /...
               (tangentParams_m * initial_current + tangentParams_e);
        
        % 处理数据并填充输出矩阵
        output_matrix(1, output_col) = Current_Value;
        output_matrix(2:513, output_col) = round(raw_data * incre);
        
        output_col = output_col + 1;  % 移动输出列指针
    end
end

% 写入单个文件
output_path = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\PLA.xlsx';
xlswrite(output_path, output_matrix);
fprintf('成功生成包含800组数据的文件：%s\n', output_path);