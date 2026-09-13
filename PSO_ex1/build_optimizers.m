function [P_main, P_para] = build_optimizers(p)
    F_B1_p = sdpvar(2, 2, 'full'); F_B2_p = sdpvar(2, 2, 'full');
    F_F1_p = sdpvar(2, 2, 'full'); F_F2_p = sdpvar(2, 2, 'full');
    F_A1_p = sdpvar(2, 2, 'full'); F_A2_p = sdpvar(2, 2, 'full');
    F_E1_p = sdpvar(2, 2, 'full'); F_E2_p = sdpvar(2, 2, 'full');
    
    parameters_in = {F_B1_p, F_B2_p, F_F1_p, F_F2_p, F_A1_p, F_A2_p, F_E1_p, F_E2_p};

    w = ones(4, 1);
    I_y = eye(2);
    O = zeros(2);
    
    ILC1 = (I_y + p.Lambda) * p.C1;
    ILC2 = (I_y + p.Lambda) * p.C2;

    hat_A = cell(1, 2); 
    
    hat_A{1} = [p.A1, O; F_A1_p - p.A1 + F_B1_p * ILC1, F_A1_p];
    hat_A{2} = [p.A2, O; F_A2_p - p.A2 + F_B2_p * ILC2, F_A2_p];

    f = cell(p.mode_number, p.divide_number + 1);  
    for i = 1 : p.mode_number
        for j = 1 : p.divide_number + 1
            f{i, j} = sdpvar(4, 1);
        end
    end
    ind = sdpvar(1);
    lar = sdpvar(1); 

    e1 = 1e-4;
    Cons_main = [];
    Cons_para = [];

    for i = 1 : p.mode_number
        for j = 1 : p.mode_number
            if i ~= j 
                Cons_main = [Cons_main, f{j, 1} - p.xi * f{i, p.divide_number + 1} - ind * w <= -e1];
                Cons_para = [Cons_para, f{j, 1} - p.xi * f{i, p.divide_number + 1} - lar * w <= -e1];
            end
        end
    end

    for i = 1 : p.mode_number
        for q = 1 : p.divide_number
            aux1 = p.divide_number * (f{i, q + 1} - f{i, q}) / p.tmin;
            
            Cons_main = [Cons_main, hat_A{i} * f{i, q} + aux1 - p.lambda * f{i, q} - ind * w <= -e1];
            Cons_main = [Cons_main, hat_A{i} * f{i, q + 1} + aux1 - p.lambda * f{i, q + 1} - ind * w <= -e1];
            
            Cons_para = [Cons_para, hat_A{i} * f{i, q} + aux1 - p.lambda * f{i, q} - lar * w <= -e1];
            Cons_para = [Cons_para, hat_A{i} * f{i, q + 1} + aux1 - p.lambda * f{i, q + 1} - lar * w <= -e1];
        end
    end

    for i = 1 : p.mode_number 
        Cons_main = [Cons_main, hat_A{i} * f{i, p.divide_number + 1} - p.lambda * f{i, p.divide_number + 1} - ind * w <= -e1];
        Cons_para = [Cons_para, hat_A{i} * f{i, p.divide_number + 1} - p.lambda * f{i, p.divide_number + 1} - lar * w <= -e1];
    end

    for i = 1 : p.mode_number
        for j = 1 : p.divide_number + 1 
            Cons_main = [Cons_main, e1 <= f{i, j} <= 2000];
            Cons_para = [Cons_para, e1 <= f{i, j} <= 2000];
        end
    end
    Cons_main = [Cons_main, ind <= 0];

    options = sdpsettings('solver', 'mosek', 'cachesolvers', 1, 'verbose', 0, ...
                          'showprogress', 0, 'mosek.MSK_IPAR_LOG', 0, ...
                          'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-7);
                          
    P_main = optimizer(Cons_main, ind, options, parameters_in, ind);
    P_para = optimizer(Cons_para, lar, options, parameters_in, lar);
end