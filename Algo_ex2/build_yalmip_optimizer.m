function P = build_yalmip_optimizer(tmin, divide_number, mode_number)
    
    %% 系统参数

    A{1} = [-3.60 1.02 0.52; 0.97 -3.40 0.82; 0.47 0.78 -3.20];
    A{2} = [-3.10 1.02 0; 0.97 -3.40 0.82; 0 0.78 -2.70];
    
    B{1} = [0.1032 0.2269 0.2159]';
    B{2} = [0.1681 0.1895 0.1376]';
    
    C{1} = [0.7165 0.8395 0.7548; 0.7266 0.7146 0.8842; 0.8848 0.7452 0.7991];
    C{2} = [0.8074 0.8001 0.7065; 0.8424 0.7494 0.7586; 0.7303 0.8909 0.8349];
    
    D{1} = [0.1798 0.1507 1.1006]';
    D{2} = [0.1442 0.1069 0.2358]';
    
    E{1} = [0.6351 0.4383 0.5213; 0.6515 0.6732 0.6919; 0.4194 0.6849 0.5771];
    E{2} = [0.6205 0.4271 0.6223; 0.5871 0.6387 0.6665; 0.6948 0.4019 0.6692];
    
    F{1} = [0.1835 0.2158 0.1677]';
    F{2} = [0.1835 0.2158 0.1677]';
    
    Lambda = diag([0.36 0.25 0.27]);
    
    F_B{1} = [0.0786 1.6002 0.0003; 0.1723 1.3878 0.0287; 0.0067 1.3721 0.0430];
    F_A{1} = [-4.6295 0.0071 0.0008; 0.0016 -4.3654 0.0044; 0.0001 0.0001 -4.2634];
    F_B{2} = [0.0175 0.2783 3.2173; 0.2456 2.6515 0.0064; 0.0003 0.1476 0.9268];
    F_A{2} = [-5.2349 0.0008 0.0019; 0.0033 -5.2349 0.0011; 0.0001 0.0005 -3.4376];
    F_F{1} = [0.0001 0.9472 0.0680; 0.0001 1.0120 0.0943; 0.2159 0.8743 0.0295];
    F_E{1} = [0.0001 0.0001 0.0001; 0.0001 0.0001 0.0001; 0.0001 0.0014 0.0001];
    F_F{2} = [0.0001 0.3889 0.5638; 0.0001 0.6742 0.2827; 0.0001 0.0667 1.0670];
    F_E{2} = [0.0002 0.0001 0.0001; 0.0001 0.0001 0.0358; 0.0002 0.0001 0.0001];
    
    I_y = eye(3);
    O = zeros(3);

    hat_A = cell(mode_number); 
    hat_B = cell(mode_number); 
    hat_C = cell(mode_number); 
    hat_D = cell(mode_number); 
    
    for i = 1 : mode_number
        hat_A{i} = [A{i} O; F_A{i} - A{i} + F_B{i} * (I_y + Lambda) * C{i} F_A{i}];
        hat_B{i} = [B{i}; F_B{i} * (I_y + Lambda) * D{i} - B{i}];
        hat_C{i} = [F_E{i} - E{i} + F_F{i} * (I_y + Lambda) * C{i} F_E{i}];
        hat_D{i} = [F_F{i} * (I_y + Lambda) * D{i} - F{i}];
    end
    
    %% 精度参数

    e1 = 1e-4;

    %% 辅助参数

    r = ones(3, 1);
    
    %% 定义 YALMIP 变量

    lambda_param = sdpvar(1, 1);
    xi_param = sdpvar(1, 1); 
    gamma_param = sdpvar(1, 1);
    
    f = cell(mode_number, divide_number + 1);  
    for i = 1 : mode_number
        for j = 1 : divide_number + 1
            f{i, j} = sdpvar(6, 1);
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
            Lvi1 = hat_A{i} * f{i, q} + aux1 - lambda_param * f{i, q} + hat_C{i}' * r;
            Constraints = [Constraints, Lvi1 <= -e1];
            Lvi2 = hat_A{i} * f{i, q + 1} + aux1 - lambda_param * f{i, q + 1} + hat_C{i}' * r;
            Constraints = [Constraints, Lvi2 <= -e1]; 
        end
    end
    
    % (2)式
    for i = 1 : mode_number 
        Lvi = hat_A{i} * f{i, divide_number + 1} - lambda_param * f{i, divide_number + 1} + hat_C{i}' * r;
        Constraints = [Constraints, Lvi <= -e1]; 
    end

    % (3)式 
    for i = 1 : mode_number 
        for q = 1 : divide_number + 1
            Lvi = hat_B{i}' * f{i, q} - gamma_param * 1 + hat_D{i}' * r;
            Constraints = [Constraints, Lvi <= -e1]; 
        end
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