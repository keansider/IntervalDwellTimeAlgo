function [U1, U2, U3] = gen_switching_signal(T_final, t_min, t_max)
%% 驻留时间信号随机生成函数
% 本函数主要用于生成两个连续子系统的切换信号
% 若需扩展至多个子系统，可在此框架基础上修改 rand() 函数的映射逻辑进行实现

%% 输入
% T_final -- 生成信号的总长度
% t_min   -- 最小驻留时间
% t_max   -- 最大驻留时间

%% 输出
% U1      -- 切换信号

% U2      -- 将每个子系统对应的驻留时间区域区分开来
%            每个子系统都有两部分，一个是[tk,tk+\tau),一个是[tk+\tau,tk+1)
%            在子系统1中第一部分信号对应1，第二部分对应2
%            在子系统2中第一部分信号对应3，第二部分对应4
%            这样做主要是为了实际时变控制器或者画V函数图时作为一个辅助识别信号

% U3      -- 时变V函数权值信号 
%            同样是为了画V函数图

%% code
    %% step1 参数校验
    if t_min >= t_max
        error('错误:请注意t_min<=t_max');
    end

    %% step2 数组预分配 
    % 预估最大切换次数，多分配 10 个空间防止越界
    max_switches = ceil(T_final / t_min) + 10;
    
    time_U1 = zeros(1, 2 * max_switches + 1);
    mode_U1 = zeros(1, 2 * max_switches + 1);
    
    time_U23 = zeros(1, 4 * max_switches);
    mode_U2  = zeros(1, 4 * max_switches);
    val_U3   = zeros(1, 4 * max_switches);

    %% step3 初始化 
    t = 0;

    % 默认第一个信号由子系统1出发
    % 但是在理论中，同时也可能由其他子系统触发
    current_mode = 1;
    
    time_U1(1) = 0;
    idx_U1  = 1;
    idx_U23 = 1;

    %% step4 主循环
    while t < T_final
        % 随机生成本次驻留时间 dt
        dt = t_min + (t_max - t_min) * rand();
        val_end = (dt - t_min) / (t_max - t_min);
        
        % --- 计算并填充 U2 和 U3 信号 (每个周期 4 个关键点) ---
        idx_range = idx_U23 : idx_U23 + 3;
        time_U23(idx_range) = [t, t + t_min, t + t_min, t + dt];
        val_U3(idx_range)   = [0, 1, 0, val_end];
        
        % --- 计算并填充 U1 信号，并分配 U2 阶段值 ---
        if current_mode == 1
            mode_U1(idx_U1 : idx_U1 + 1) = [1, 1];
            mode_U2(idx_range) = [1, 1, 2, 2];
            next_mode = 2;
        else
            mode_U1(idx_U1 : idx_U1 + 1) = [2, 2];
            mode_U2(idx_range) = [3, 3, 4, 4];
            next_mode = 1;
        end
        
        % 更新时间
        t = t + dt;
        
        idx_U1 = idx_U1 + 2;
        time_U1(idx_U1 - 1 : idx_U1) = [t, t];
        
        % 更新索引与下一阶段模式
        idx_U23 = idx_U23 + 4;
        current_mode = next_mode;
    end

    %% step5 截断多余数组与收尾处理
    % 闭合最后一个状态的模式
    mode_U1(idx_U1) = current_mode;
    
    time_U1 = time_U1(1 : idx_U1);
    mode_U1 = mode_U1(1 : idx_U1);
    
    time_U23 = time_U23(1 : idx_U23 - 1);
    mode_U2  = mode_U2(1 : idx_U23 - 1);
    val_U3   = val_U3(1 : idx_U23 - 1);

    %% step6 输出结果
    U1 = [time_U1; mode_U1]';
    U2 = [time_U23; mode_U2]';
    U3 = [time_U23; val_U3]';
end