% KL散度计算函数
function kl = kldiv(p, q)
    epsilon = 1e-10; % 避免零值
    p = (p + epsilon) / (sum(p) + epsilon*numel(p));
    q = (q + epsilon) / (sum(q) + epsilon*numel(q));
    kl = sum(p .* log(p ./ q));
end