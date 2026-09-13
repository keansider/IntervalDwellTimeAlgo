function [P_main, P_para] = build_optimizers(p)
    disp('正在预编译 YALMIP 模型 (Main 和 Para)... 这只需执行一次，请稍候...');
    
    %% 1. 将待优化的滤波器矩阵定义为 YALMIP 参数 (Inputs)
    F_B1_p = sdpvar(3, 3, 'full'); F_B2_p = sdpvar(3, 3, 'full');
    F_F1_p = sdpvar(3, 3, 'full'); F_F2_p = sdpvar(3, 3, 'full');
    F_A1_p = sdpvar(3, 3, 'full'); F_A2_p = sdpvar(3, 3, 'full');
    F_E1_p = sdpvar(3, 3, 'full'); F_E2_p = sdpvar(3, 3, 'full');
    
    % 注意：传入和传出的顺序必须与 evaluate_particle_fitness 中的 inputs 顺序严格对应
    parameters_in = {F_B1_p, F_B2_p, F_F1_p, F_F2_p, F_A1_p, F_A2_p, F_E1_p, F_E2_p};

    %% 2. 定义辅助矩阵与系统常量
    w = ones(6, 1);
    r = ones(3, 1);
    I_y = eye(3);
    O = zeros(3);
    
    % 提取公共计算以加速编译
    ILC1 = (I_y + p.Lambda) * p.C1;
    ILD1 = (I_y + p.Lambda) * p.D1;
    ILC2 = (I_y + p.Lambda) * p.C2;
    ILD2 = (I_y + p.Lambda) * p.D2;

    %% 3. 构建包含参数的系统矩阵 hat_A, hat_B, hat_C, hat_D
    hat_A = cell(1, 2); hat_B = cell(1, 2); hat_C = cell(1, 2); hat_D = cell(1, 2);
    
    hat_A{1} = [p.A1, O; F_A1_p - p.A1 + F_B1_p * ILC1, F_A1_p];
    hat_A{2} = [p.A2, O; F_A2_p - p.A2 + F_B2_p * ILC2, F_A2_p];

    hat_B{1} = [p.B1; F_B1_p * ILD1 - p.B1];
    hat_B{2} = [p.B2; F_B2_p * ILD2 - p.B2];

    hat_C{1} = [F_E1_p - p.E1 + F_F1_p * ILC1, F_E1_p];
    hat_C{2} = [F_E2_p - p.E2 + F_F2_p * ILC2, F_E2_p];

    hat_D{1} = F_F1_p * ILD1 - p.F1;
    hat_D{2} = F_F2_p * ILD2 - p.F2;

    %% 4. 定义决策变量 (Decision Variables)
    f = cell(p.mode_number, p.divide_number + 1);  
    for i = 1 : p.mode_number
        for j = 1 : p.divide_number + 1
            f{i, j} = sdpvar(6, 1);
        end
    end
    gamma = sdpvar(1);
    lar = sdpvar(1); 

    %% 5. 同步构建主模型(Main)和罚函数模型(Para)的约束
    e1 = 1e-4;
    Cons_main = [];
    Cons_para = [];

    % (1)式
    for i = 1 : p.mode_number
        for j = 1 : p.mode_number
            if i ~= j 
                Cons_main = [Cons_main, f{j, 1} - p.xi * f{i, p.divide_number + 1} <= -e1];
                Cons_para = [Cons_para, f{j, 1} - p.xi * f{i, p.divide_number + 1} <= -e1];
            end
        end
    end

    % (2).1式 与 (2).2式
    for i = 1 : p.mode_number
        for q = 1 : p.divide_number
            aux1 = p.divide_number * (f{i, q + 1} - f{i, q}) / p.tmin;
            
            Cons_main = [Cons_main, hat_A{i}' * f{i, q} + aux1 - p.lambda * f{i, q} + hat_C{i}' * r <= -e1];
            Cons_main = [Cons_main, hat_A{i}' * f{i, q + 1} + aux1 - p.lambda * f{i, q + 1} + hat_C{i}' * r <= -e1];
            
            Cons_para = [Cons_para, hat_A{i}' * f{i, q} + aux1 - p.lambda * f{i, q} + hat_C{i}' * r - lar * w <= -e1];
            Cons_para = [Cons_para, hat_A{i}' * f{i, q + 1} + aux1 - p.lambda * f{i, q + 1} + hat_C{i}' * r - lar * w <= -e1];
        end
    end

    % (2)式 末端
    for i = 1 : p.mode_number 
        Cons_main = [Cons_main, hat_A{i}' * f{i, p.divide_number + 1} - p.lambda * f{i, p.divide_number + 1} + hat_C{i}' * r <= -e1];
        Cons_para = [Cons_para, hat_A{i}' * f{i, p.divide_number + 1} - p.lambda * f{i, p.divide_number + 1} + hat_C{i}' * r - lar * w <= -e1];
    end

    % (3)式
    for i = 1 : p.mode_number 
        for q = 1 : p.divide_number + 1
            Cons_main = [Cons_main, hat_B{i}' * f{i, q} - gamma + hat_D{i}' * r <= -e1];
            Cons_para = [Cons_para, hat_B{i}' * f{i, q} - gamma + hat_D{i}' * r - lar * r <= -e1]; 
        end
    end

    % 正约束 与 辅助约束
    for i = 1 : p.mode_number
        for j = 1 : p.divide_number + 1 
            Cons_main = [Cons_main, e1 <= f{i, j} <= 1000];
            Cons_para = [Cons_para, e1 <= f{i, j} <= 1000];
        end
    end
    Cons_main = [Cons_main, gamma >= e1];
    Cons_para = [Cons_para, gamma >= e1];

    %% 6. 编译优化器对象
    options = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0, ...
                          'showprogress', 0, 'mosek.MSK_IPAR_LOG', 0, ...
                          'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-7);
                          
    outputs_main = gamma;  % 【修改】去掉大括号，强制要求按普通数值(double)输出
    outputs_para = lar;    % 【修改】去掉大括号
    
    P_main = optimizer(Cons_main, gamma, options, parameters_in, outputs_main);
    P_para = optimizer(Cons_para, lar, options, parameters_in, outputs_para);
    
    disp('模型预编译完成！PSO 正式开始迭代...');
end