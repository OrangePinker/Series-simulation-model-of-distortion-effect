% SNR计算函数
function snr = calc_snr(signal)
    window_size = 5;
    smoothed = movmean(signal, window_size);
    noise = signal - smoothed;
    noise_std = std(noise);
    % 防止分母为零
    if noise_std < eps
        noise_std = eps;  % 或根据需求设置阈值
    end
    snr = 20*log10(mean(signal)/noise_std);

end