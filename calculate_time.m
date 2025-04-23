function t = calculate_time(Q0, Q1, U, d)
    % 设置默认参数值
    if nargin < 1
        Q0 = 8.8539e-6;  % Q0 (C/m^3)
    end
    if nargin < 2
        Q1 = 1e-6;       % Q1 (C/m^3)
    end
    if nargin < 3
        U = 500;         % 电压 (V)
    end
    if nargin < 4
        d = 0.01;        % 运动距离 (m)
    end
    
    % 定义常数
    Lambda = 0.05;      % Λ (m)
    epsilon0 = 8.854e-12; % 真空介电常数 (F/m)
    q = -1.6e-19;       % 负电荷量 (C)
    m = 9.11e-31;       % 电子质量 (kg)
    
    % 计算电场强度
    E = -U / d;          % 电场强度 (V/m), 方向从正极指向负极
    
    % 定义电势的一阶导数 dV/dz
    dVdz = @(z) (Q0 * Lambda * exp(-z/Lambda)) / epsilon0 - (Q1 * z) / epsilon0;
    
    % 定义电荷分布产生的力 F_distribution(z) = q * dV/dz（反向）
    F_distribution = @(z) q * dVdz(z);
    
    % 定义匀强电场产生的力 F_external = q * E（方向不变）
    F_external = q * E;
    
    % 定义总力 F_total(z) = F_distribution(z) + F_external
    F_total = @(z) F_distribution(z) + F_external;
    
    % 定义运动方程 d^2z/dt^2 = F_total(z)/m
    d2zdt2 = @(t, z) [z(2); F_total(z(1))/m]; % z(1) = z, z(2) = dz/dt
    
    % 初始条件
    z0 = 0;             % 初始位置（从负极出发）
    v0 = 0;             % 初始速度
    initial_conditions = [z0; v0];
    
    % 时间范围（设置为足够大的范围但不超过必要时间）
    tspan = [0 1e-3];   % 时间范围 (s)
    
    % 使用 ode45 求解运动方程
    options = odeset('Events', @(t, z) event_function(t, z, d), 'RelTol', 1e-6, 'AbsTol', 1e-9);
    [t, z, te, ze] = ode45(d2zdt2, tspan, initial_conditions, options);
    
    % 提取运动时间
    if ~isempty(te)
        t = te;  % 运动时间
        fprintf('运动到距离 d = %.4f m 所用时间 t = %.6e s\n', d, t);
    else
        t = NaN; % 未能在给定时间内达到距离 d
        fprintf('未能在给定时间内达到距离 d = %.4f m\n', d);
    end
    

    % 事件函数：当 z 达到 d 时停止计算
    function [value, isterminal, direction] = event_function(~, z, d)
        value = z(1) - d; % 当 z(1) = d 时触发事件
        isterminal = 1;   % 停止计算
        direction = 1;    % 仅当 z(1) 增加时触发
    end
end