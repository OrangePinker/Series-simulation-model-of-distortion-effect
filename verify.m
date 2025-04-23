%% 数据读取
real_data = readtable('v1.xlsx', 'Range', 'A2:T513');  % 根据实际列范围调整BO
sim_data = readtable('a1.xlsx', 'Range', 'A2:T513');
real_data = table2array(real_data);
sim_data = table2array(sim_data);

%% 条件1：全局KL散度 & KS检验

% 分布对比
[real_counts, edges] = histcounts(real_data(:), 100);
sim_counts = histcounts(sim_data(:), edges);
kl_value = kldiv(real_counts, sim_counts);

% KS检验
[~, p_ks] = kstest2(real_data(:), sim_data(:));
condition1 = (kl_value < 0.1) && (p_ks > 0.05);

% 输出条件1的结果
disp(['条件1: KL散度 = ', num2str(kl_value), ', KS检验p值 = ', num2str(p_ks)]);

%% 条件2：特征峰位置 & SNR差异
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
condition2 = (peak_shift <= 1) && (snr_diff < 0.1);

% 输出条件2的结果
disp(['条件2: 主峰位移 = ', num2str(peak_shift), ', SNR差异 = ', num2str(snr_diff)]);


%% 最终判定
if condition1 && condition2 
    disp('仿真数据验证通过，可用于替代实际数据');
else
    failed_cond = find([condition1, condition2] == 0);
    disp(['验证失败条件：', num2str(failed_cond)]);
end