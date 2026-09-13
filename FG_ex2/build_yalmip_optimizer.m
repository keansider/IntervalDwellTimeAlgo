function P = build_yalmip_optimizer(tmin, divide_number, mode_number)
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
    F{2} = [0.1457 0.1087 0.2149]';
    Lambda = diag([0.26 0.15 0.17]);
    F_B{1} = [1.4178 0.0193 0.0183; 0.9018 0.8022 0.0131; 0.2699 1.2881 0.0200];
    F_A{1} = [-4.3655 0.3882 0.0199; 0.0405 -4.1208 0.0058; 0.0035 0.0021 -4.2692];
    F_B{2} = [3.1045 0.6203 0.0013; 1.6001 1.2605 0.0652; 0.2473 0.6981 0.2858];
    F_A{2} = [-5.2349 0.0001 2.1672; 0.0174 -5.1705 0.0030; 0.0055 0.0001 -3.4341];
    F_F{1} = [0.7072 0.3659 0.0467; 0.0212 0.9892 0.0952; 0.0551 1.0297 0.0315];
    F_E{1} = [0.0003 0.0001 0.0004; 0.0001 0.0017 0.0001; 0.0004 0.0070 0.0014];
    F_F{2} = [0.0057 0.3999 0.5574; 0.0003 0.3264 0.5825; 0.0022 0.0607 1.0726];
    F_E{2} = [0.0001 0.0001 0.0001; 0.0016 0.0003 0.0529; 0.0028 0.0012 0.0003];
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
    e1 = 1e-4;
    r = ones(3, 1);
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
    for i = 1 : mode_number
        for j = 1 : mode_number
            if i ~= j 
                Lvi = f{j, 1} - xi_param * f{i, divide_number + 1};
                Constraints = [Constraints, Lvi <= -e1];
            end
        end
    end
    for i = 1 : mode_number
        for q = 1 : divide_number
            aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
            Lvi1 = hat_A{i}' * f{i, q} + aux1 - lambda_param * f{i, q} + hat_C{i}' * r;
            Constraints = [Constraints, Lvi1 <= -e1];
            Lvi2 = hat_A{i}' * f{i, q + 1} + aux1 - lambda_param * f{i, q + 1} + hat_C{i}' * r;
            Constraints = [Constraints, Lvi2 <= -e1]; 
        end
    end
    for i = 1 : mode_number 
        Lvi = hat_A{i}' * f{i, divide_number + 1} - lambda_param * f{i, divide_number + 1} + hat_C{i}' * r;
        Constraints = [Constraints, Lvi <= -e1]; 
    end
    for i = 1 : mode_number 
        for q = 1 : divide_number + 1
            Lvi = hat_B{i}' * f{i, q} - gamma_param * 1 + hat_D{i}' * r;
            Constraints = [Constraints, Lvi <= -e1]; 
        end
    end 
    for i = 1 : mode_number
        for j = 1 : divide_number + 1 
            Constraints = [Constraints, f{i, j} >= e1]; 
        end
    end
    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'mosek.MSK_IPAR_LOG', 0, 'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-7, 'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-7, 'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-7);
    P = optimizer(Constraints, gamma_param, options, [lambda_param; xi_param], gamma_param);
end