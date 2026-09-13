function P = build_yalmip_optimizer(tmin, divide_number, mode_number)
    
    %% 系统参数

    A{1} = [-2.487 0.231; 0.322 0.402];
    A{2} = [0.366 0.311; 0.392 -2.137];
    B{1} = [0.143 0.223]';
    B{2} = [0.147 0.303]';
    C{1} = [0.758 0.826; 0.822 0.874];
    C{2} = [0.898 0.842; 0.762 0.827];
    D{1} = [0.145 0.269]';
    D{2} = [0.125 0.249]';
    E{1} = [0.658 0.682; 0.694 0.427];
    E{2} = [0.648 0.462; 0.623 0.617];
    F{1} = [0.225 0.261]';
    F{2} = [0.185 0.279]';
    Lambda = diag([0.1, 0.2]);
    
    F_B{1} = [3.0160 0.0029; 0.5439 0.0256];
    F_B{2} = [1.4346 0.0303; 2.7971 0.0121];
    F_A{1} = [-4.5377 0.0035; 0.0270 -0.0168];
    F_A{2} = [-0.8108 0.0059; 0.0103 -4.2193];
    F_F{1} = [3.8209 3.8209; 3.7631 3.8209];
    F_F{2} = [1.3228 3.0085; 3.8209 3.8209];
    F_E{1} = [2.5763 3.8209; 3.8209 1.1209];
    F_E{2} = [3.8209 3.8209; 0.8349 3.8209];
    
    I_y = eye(2);
    O = zeros(2);
    
    for i = 1 : mode_number
        hat_A{i} = [A{i} O; F_A{i} - A{i} + F_B{i} * (I_y + Lambda) * C{i} F_A{i}];
    end
    
    %% 精度参数

    e1 = 1e-4;
    
    %% 定义 YALMIP 变量

    lambda_param = sdpvar(1, 1);
    xi_param = sdpvar(1, 1);
    
    f = cell(mode_number, divide_number + 1);  
    for i = 1 : mode_number
        for j = 1 : divide_number + 1
            f{i, j} = sdpvar(4, 1);
        end
    end
    
    Constraints = [];
    
    %% 构造约束 

    % (1)式
    for i = 1 : mode_number
        for j = 1 : mode_number
            if i ~= j 
                Lvi = f{j, 1} - xi_param * f{i, divide_number + 1};
                Constraints = [Constraints, Lvi <= -e1];
            end
        end
    end
    
    % (2).1式 & (2).2式
    for i = 1 : mode_number
        for q = 1 : divide_number
            aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
            Lvi1 = hat_A{i} * f{i, q} + aux1 - lambda_param * f{i, q};
            Constraints = [Constraints, Lvi1 <= -e1];
            Lvi2 = hat_A{i} * f{i, q + 1} + aux1 - lambda_param * f{i, q + 1};
            Constraints = [Constraints, Lvi2 <= -e1]; 
        end
    end
    
    % (2)式
    for i = 1 : mode_number 
        Lvi = hat_A{i} * f{i, divide_number + 1} - lambda_param * f{i, divide_number + 1};
        Constraints = [Constraints, Lvi <= -e1]; 
    end
    
    % 正约束与辅助约束
    for i = 1 : mode_number
        for j = 1 : divide_number + 1 
            Constraints = [Constraints, f{i, j} >= e1, f{i, j} <= 1000]; 
        end
    end
    
    %% 求解器设置
    options = sdpsettings('solver', 'mosek', ...
                          'verbose', 0, ...
                          'mosek.MSK_IPAR_LOG', 0, ... 
                          'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-7,...
                          'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-7,...
                          'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-7);
                          
    %% 【关键】生成 optimizer 对象
    % 输入：[lambda_param; xi_param] (必须是列向量)
    % 输出：随便指定一个，这里输出 f{1,1}(1) 仅为占位，因为我们只关心是否可行 (errorcode)
    P = optimizer(Constraints, 0, options, [lambda_param; xi_param], f{1,1}(1));
end