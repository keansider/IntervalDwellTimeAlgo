function [U1, U2] = gen_periodic_switching_signal(T_final, t_min)
%% 周期驻留时间信号生成函数
% 本函数用于生成两个连续子系统的周期切换信号及对应的斜率信号
%% 输入
% T_final -- 生成信号的总长度
% t_min   -- 固定的驻留时间（周期时间）
%% 输出
% U1      -- 切换信号 (方波：在子系统 1 和 2 之间交替)
% U2      -- 斜率信号 (锯齿波：在每个 t_min 时间间隔内从 0 线性增加至 1)
%% code
    %% step1 参数校验
    if t_min <= 0
        error('错误:请注意固定的驻留时间 t_min 必须大于0');
    end
    
    %% step2 数组预分配 
    % 预估最大切换次数，多分配 2 个空间防止越界
    max_switches = ceil(T_final / t_min) + 2;
    
    % 每个周期需要2个关键点来绘制垂直跳变
    time_U = zeros(1, 2 * max_switches);
    mode_U = zeros(1, 2 * max_switches);
    val_U  = zeros(1, 2 * max_switches);
    
    %% step3 初始化 
    t = 0;
    % 默认第一个信号由子系统1出发
    current_mode = 1; 
    idx = 1;
    
    %% step4 主循环
    while t < T_final
        % 下一个切换点的时间
        t_next = t + t_min;
        
        % --- 计算并填充 U1 和 U2 信号 (每个周期 2 个关键点) ---
        % U1: 保持 current_mode 直到 t_next
        % U2: 从 0 线性增加到 1
        time_U(idx : idx + 1) = [t, t_next];
        mode_U(idx : idx + 1) = [current_mode, current_mode];
        val_U(idx : idx + 1)  = [0, 1];
        
        % 更新时间
        t = t_next;
        
        % 交替切换子系统状态
        if current_mode == 1
            current_mode = 2;
        else
            current_mode = 1;
        end
        
        % 更新索引
        idx = idx + 2;
    end
    
    %% step5 截断多余数组与收尾处理
    time_U = time_U(1 : idx - 1);
    mode_U = mode_U(1 : idx - 1);
    val_U  = val_U(1 : idx - 1);
    
    %% step6 输出结果
    U1 = [time_U; mode_U]';
    U2 = [time_U; val_U]';
end