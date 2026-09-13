function [pop, gbest_pos, gbest_fitness, gbest_penalty, gbest_is_valid] = init_population(p)
    %% 初始化结构体与变量 
    pop = struct('pos', {}, 'vel', {}, ...
                 'is_valid', {}, 'fitness', {}, 'penalty', {}, ...
                 'pbest_pos', {}, 'pbest_is_valid', {}, ...
                 'pbest_fitness', {}, 'pbest_penalty', {}); 
             
    gbest_pos = [];
    gbest_is_valid = false;
    gbest_fitness = inf;  
    gbest_penalty = inf; 
    HUGE_PENALTY = 1e8; 
    
    fprintf('开始极速初始化种群，总数: %d\n', p.pop_size);
    
    % --- 【提速核心】：提取公共常数项矩阵，避免循环内重复计算 ---
    I_y = eye(3);
    ILC1 = (I_y - p.Lambda) * p.C1;
    ILD1 = (I_y - p.Lambda) * p.D1;
    ILC2 = (I_y - p.Lambda) * p.C2;
    ILD2 = (I_y - p.Lambda) * p.D2;
    mask_off_diag = ~I_y; % 非对角线掩码
    % -------------------------------------------------------------

    for i = 1 : p.pop_size 
        %% 1. 初始化 A1 B1 (边界推导法)
        is_AB1_valid = false;
        while ~is_AB1_valid
            % 第 1 步：生成 B1 (必须 >= 0 且满足 m2 >= 0)
            F_B1 = p.x_max * rand(p.dims.F_B1);
            if any(any( F_B1 * ILD1 - p.B1 < 0 ))
                continue; % 不满足 m2 则重刷 B1
            end
            
            % 第 2 步：计算 A1 的下界 (保证 m1 >= 0)
            LB_A = p.A1 - F_B1 * ILC1;
            LB_A(mask_off_diag) = max(0, LB_A(mask_off_diag)); % 非对角线必须 >= 0
            LB_A(~mask_off_diag) = max(p.x_min, LB_A(~mask_off_diag)); % 对角线受 x_min 限制
            
            if any(LB_A(:) > p.x_max)
                continue; % 如果推导出的下界超出了全局上界 x_max，说明这个 B1 是死胡同，重刷 B1
            end
            
            % 第 3 步：在安全区间内直接生成 A1 (命中率 100%)
            F_A1 = LB_A + (p.x_max - LB_A) .* rand(p.dims.F_A1);
            
            pop(i).pos.F_B1 = F_B1; pop(i).vel.F_B1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_B1);
            pop(i).pos.F_A1 = F_A1; pop(i).vel.F_A1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_A1);
            is_AB1_valid = validate_AB(pop(i).pos, p, 1); % 保底验证，理论上必然通过
        end
        
        %% 2. 初始化 A2 B2 (边界推导法)
        is_AB2_valid = false;
        while ~is_AB2_valid
            F_B2 = p.x_max * rand(p.dims.F_B2);
            if any(any( F_B2 * ILD2 - p.B2 < 0 ))
                continue;
            end
            
            LB_A = p.A2 - F_B2 * ILC2;
            LB_A(mask_off_diag) = max(0, LB_A(mask_off_diag)); 
            LB_A(~mask_off_diag) = max(p.x_min, LB_A(~mask_off_diag)); 
            
            if any(LB_A(:) > p.x_max)
                continue;
            end
            
            F_A2 = LB_A + (p.x_max - LB_A) .* rand(p.dims.F_A2);
            
            pop(i).pos.F_B2 = F_B2; pop(i).vel.F_B2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_B2);
            pop(i).pos.F_A2 = F_A2; pop(i).vel.F_A2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_A2);
            is_AB2_valid = validate_AB(pop(i).pos, p, 2); 
        end
        
        %% 3. 初始化 E1 F1 (边界推导法)
        is_EF1_valid = false;
        while ~is_EF1_valid
            F_F1 = p.x_max * rand(p.dims.F_F1);
            if any(any( F_F1 * ILD1 - p.F1 < 0 ))
                continue;
            end
            
            LB_E = p.E1 - F_F1 * ILC1;
            LB_E = max(0, LB_E); % E1 的所有元素必须 >= 0
            
            if any(LB_E(:) > p.x_max)
                continue;
            end
            
            F_E1 = LB_E + (p.x_max - LB_E) .* rand(p.dims.F_E1);
            
            pop(i).pos.F_F1 = F_F1; pop(i).vel.F_F1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_F1);
            pop(i).pos.F_E1 = F_E1; pop(i).vel.F_E1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_E1);
            is_EF1_valid = validate_EF(pop(i).pos, p, 1);
        end
        
        %% 4. 初始化 E2 F2 (边界推导法)
        is_EF2_valid = false;
        while ~is_EF2_valid
            F_F2 = p.x_max * rand(p.dims.F_F2);
            if any(any( F_F2 * ILD2 - p.F2 < 0 ))
                continue;
            end
            
            LB_E = p.E2 - F_F2 * ILC2;
            LB_E = max(0, LB_E); % E2 的所有元素必须 >= 0
            
            if any(LB_E(:) > p.x_max)
                continue;
            end
            
            F_E2 = LB_E + (p.x_max - LB_E) .* rand(p.dims.F_E2);
            
            pop(i).pos.F_F2 = F_F2; pop(i).vel.F_F2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_F2);
            pop(i).pos.F_E2 = F_E2; pop(i).vel.F_E2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_E2);
            is_EF2_valid = validate_EF(pop(i).pos, p, 2);
        end
        
        %% fitness 评价 
        pop(i).pbest_is_valid = false; 
        pop(i) = evaluate_particle_fitness(pop(i), p, HUGE_PENALTY);
        pop(i) = update_individual_best(pop(i));
        [gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
            update_global_best(pop(i), gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty);
            
        %% [修改处] 实时打印输出，使用标准逐行打印
        percentDone = 100 * i / p.pop_size;
        fprintf('当前进度: %d / %d (%.1f%%)\n', i, p.pop_size, percentDone);
    end
    fprintf('种群初始化完成！\n');
end