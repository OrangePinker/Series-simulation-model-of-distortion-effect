%%
% 串联模型
clc; clear;

% 记录运行开始时间
tic;

% 初始化进度条
processed_counts = 0;
percent_step = 1;
last_update_percent = 0;
h = waitbar(0,'Processing...');

% 提取增量和衰减参数（来源于总光子数）
Total_filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\averaged_result.xlsx';
x0 = 200;
[fitParams, tangentParams] = fitAndTangent(Total_filename, x0);
fitParams_a = fitParams(1);
fitParams_b = fitParams(2);
fitParams_c = fitParams(3);
fitParams_d = fitParams(4);
tangentParams_m = tangentParams(1);
tangentParams_e = tangentParams(2);

%提取分布参数（来源于未畸变光谱）
filename = 'C:\Users\Lenovo\Desktop\毕设\光子计数探测器谱线失真机理及校正方法研究\正向仿真\raw200.xlsx';
[params1, params2] = fitTwoWeibull(filename);
scale_param = params1;
shape_param = params2;

% 导入Excel文件
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

%-------------------------------增量计算------------------------------------

% 计算增量
initial_current = 200;
Current_Value = 360;  %电流值，实验时唯一变量
incre = (tangentParams_m * Current_Value + tangentParams_e) / (tangentParams_m * initial_current + tangentParams_e);

% 将导入的数据依次与incre相乘
scaled_data = zeros(1, num_bins);  % 创建存储结果的数组
for i = 1:num_bins
    scaled_data(i) = bins{i} * incre;  % 乘以incre
end

scaled_data = round(scaled_data);
%-------------------------------极化效应------------------------------------

% 极化参数设置
Upper_Limit = 600;
Lower_Limit = 200;
Polarization_Q0 = (8.84e-6 - 8.6e-6)* (Current_Value - Lower_Limit) / (Upper_Limit - Lower_Limit)+8.6e-6;
t = calculate_time(Polarization_Q0, 1e-6, 500, 0.01);
K = 2.5e6;%3e6
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
        pending_data = scaled_data(i);  % 存入新的original_data(i)
        last_non_zero_index = new_index;
    else
        pola_data(new_index) = pola_data(new_index) + scaled_data(i);

    end
end

    new_index = round(new_index);
% 极化数组大小
pola_data = round(pola_data);
sum_pola = sum(pola_data);

%-------------------------------脉冲堆积------------------------------------

% 新数组初始化
pulse_array = zeros(1, 512);

% 定义每次取出数量的概率分布
% 例如，90%的概率取1，10%的概率取2
Probability = -0.029+(fitParams_a * Current_Value^3 + fitParams_b * Current_Value^2 + fitParams_c * Current_Value + fitParams_d)/(tangentParams_m * Current_Value + tangentParams_e)
prob_take_1 = Probability;
prob_take_2 = Probability+(1 - Probability)*0.81;%0.81
prob_take_3 = 1.00;

% 循环直到scaled_data中的所有数被取完
while any(pola_data > 0)

    % 更新processed_counts
    zeroCount = sum(pola_data == 0);

    % 限制单次运行时长
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

    if bin_index < 7 || bin_index > 512 || pola_data(bin_index) < 1
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

        if bin_index2 < 1 || bin_index2 > 512 || pola_data(bin_index2) < 1
            continue;  
        end

        break;

    end

    % 检查取出的数量是否大于1
    if take_num == 1
        % 取出一个数
        if pola_data(bin_index) > 0
            pulse_array(bin_index) = pulse_array(bin_index) + 1;
            pola_data(bin_index) = pola_data(bin_index) - 1;
        end

    elseif take_num == 1.5
        compare_bin = max(bin_index2,bin_index);
        if compare_bin < 7
            continue;  
        end
        if pola_data(compare_bin) > 0
            pulse_array(bin_index) = pulse_array(bin_index) + 1;
            pola_data(bin_index) = pola_data(bin_index) - 1;
            pola_data(bin_index2) = pola_data(bin_index2) - 1;
        end
        
    else
        % 随机选择take_num个数
        selected_indices = bin_index2 + bin_index;
        if selected_indices < 7
            continue;  
        end
        
        if selected_indices > 512
            pola_data(bin_index2) = pola_data(bin_index2) - 1;
            pola_data(bin_index) = pola_data(bin_index) - 1;
            continue; 
        end
        if selected_indices <= 512
            % 检查这些bin是否有足够的数量
            sufficient = true;
            for i = 1:take_num
                if pola_data(bin_index2) < 1 || pola_data(bin_index) < 1
                    sufficient = false;
                    break;
                end
            end
            if sufficient
                pulse_array(selected_indices) = pulse_array(selected_indices) + 1;
                pola_data(bin_index2) = pola_data(bin_index2) - 1;
                pola_data(bin_index) = pola_data(bin_index) - 1;
                
            end
        end
    end
end

% 关闭进度条
close(h);
% 输出运行时间
fprintf('Total runtime: %f seconds.\n', toc);

% 归一化数据
original_data_norm = (original_data - min(original_data)) / (max(original_data) - min(original_data));
pulse_array_norm = (pulse_array - min(pulse_array)) / (max(pulse_array) - min(pulse_array));

% 计算并输出数组的和
sum_original = sum(original_data);
sum_scaled = sum(scaled_data);
sum_pulse = sum(pulse_array);

disp(['原始数据的和: ', num2str(sum_original)]);
disp(['增量数据的和: ', num2str(sum_scaled)]);
disp(['极化数据的和: ', num2str(sum_pola)]);
disp(['堆积数据的和: ', num2str(sum_pulse)]);

% 绘制图形
figure;
plot(1:512, original_data_norm, 'r-', 'DisplayName', 'Original Data','LineWidth', 1.5);
hold on;
plot(1:512, pulse_array_norm, 'b-', 'DisplayName', 'Attenuated Data','LineWidth', 1.5);

legend show;
xlabel('Index');
ylabel('Value');
title('Comparison of Original and Attenuated Data');
grid on;

% 确保 pulse_array 是列向量（512x1）
pulse_array = pulse_array(:); % 若为行向量则自动转置为列

% 写入Excel，从A2开始填充
% writematrix(pulse_array, 'a1.xlsx', 'Range', 'A2');
sim_data = pulse_array;

%-------------------------------验证部分------------------------------------



%% 数据读取
real_data = readtable('v1.xlsx', 'Range', 'A2:T513');  % 根据实际列范围调整BO
% sim_data = readtable('a1.xlsx', 'Range', 'A2:T513');
real_data = table2array(real_data);
% sim_data = table2array(sim_data);

% 绘制图形
figure;
plot(1:512, real_data, 'r-', 'DisplayName', 'Original Data','LineWidth', 1.5);
hold on;
plot(1:512, sim_data, 'b-', 'DisplayName', 'Attenuated Data','LineWidth', 1.5);

legend show;
xlabel('Index');
ylabel('Value');
title('Comparison of Original and Attenuated Data');
grid on;

%% 条件1：全局KL散度 & KS检验

% 分布对比
[real_counts, edges] = histcounts(real_data(:), 100);
sim_counts = histcounts(sim_data(:), edges);
kl_value = kldiv(real_counts, sim_counts);

% KS检验
[~, p_ks] = kstest2(real_data(:), sim_data(:));
condition1 = (kl_value < 0.5) && (p_ks > 0.05);

% 输出条件1的结果
disp(['条件1: KL散度 = ', num2str(kl_value), ', KS检验p值 = ', num2str(p_ks)]);

%% 条件2：特征峰位置 
% 平均光谱分析
mean_real = mean(real_data, 1);
mean_sim = mean(sim_data, 1);

% 主峰定位
[~, real_loc] = max(mean_real);
[~, sim_loc] = max(mean_sim);
peak_shift = abs(real_loc - sim_loc);

snr_real = calc_snr(mean_real);
snr_sim = calc_snr(mean_sim);
snr_diff = abs(snr_real - snr_sim) / snr_real;
condition2 = (peak_shift <= 1) ;

% 输出条件2的结果
disp(['条件2: 主峰位移 = ', num2str(peak_shift), ]);


%% 最终判定
if condition1 && condition2 
    disp('仿真数据验证通过，可用于替代实际数据');
else
    failed_cond = find([condition1, condition2] == 0);
    disp(['验证失败条件：', num2str(failed_cond)]);
end