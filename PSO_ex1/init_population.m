function [pop, gbest_pos, gbest_fitness, gbest_penalty, gbest_is_valid] = init_population(p)
    pop = struct('pos', {}, 'vel', {}, ...
                 'is_valid', {}, 'fitness', {}, 'penalty', {}, ...
                 'pbest_pos', {}, 'pbest_is_valid', {}, ...
                 'pbest_fitness', {}, 'pbest_penalty', {}); 
             
    gbest_pos = [];
    gbest_is_valid = false;
    gbest_fitness = inf;  
    gbest_penalty = inf; 
    
    HUGE_PENALTY = 1e8; 
    
    I_y = eye(2);
    ILC1 = (I_y - p.Lambda) * p.C1;
    ILD1 = (I_y - p.Lambda) * p.D1;
    ILC2 = (I_y - p.Lambda) * p.C2;
    ILD2 = (I_y - p.Lambda) * p.D2;
    mask_off_diag = ~I_y;

    for i = 1 : p.pop_size 
        is_AB1_valid = false;
        while ~is_AB1_valid
            F_B1 = p.x_max * rand(p.dims.F_B1);
            if any(any( F_B1 * ILD1 - p.B1 < 0 ))
                continue; 
            end
            
            LB_A = p.A1 - F_B1 * ILC1;
            LB_A(mask_off_diag) = max(0, LB_A(mask_off_diag)); 
            LB_A(~mask_off_diag) = max(p.x_min, LB_A(~mask_off_diag)); 
            
            if any(LB_A(:) > p.x_max)
                continue; 
            end
            
            F_A1 = LB_A + (p.x_max - LB_A) .* rand(p.dims.F_A1);
            
            pop(i).pos.F_B1 = F_B1; pop(i).vel.F_B1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_B1);
            pop(i).pos.F_A1 = F_A1; pop(i).vel.F_A1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_A1);
            is_AB1_valid = validate_AB(pop(i).pos, p, 1); 
        end
        
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
        
        is_EF1_valid = false;
        while ~is_EF1_valid
            F_F1 = p.x_max * rand(p.dims.F_F1);
            if any(any( F_F1 * ILD1 - p.F1 < 0 ))
                continue;
            end
            
            LB_E = p.E1 - F_F1 * ILC1;
            LB_E = max(0, LB_E); 
            
            if any(LB_E(:) > p.x_max)
                continue;
            end
            
            F_E1 = LB_E + (p.x_max - LB_E) .* rand(p.dims.F_E1);
            
            pop(i).pos.F_F1 = F_F1; pop(i).vel.F_F1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_F1);
            pop(i).pos.F_E1 = F_E1; pop(i).vel.F_E1 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_E1);
            is_EF1_valid = validate_EF(pop(i).pos, p, 1);
        end
        
        is_EF2_valid = false;
        while ~is_EF2_valid
            F_F2 = p.x_max * rand(p.dims.F_F2);
            if any(any( F_F2 * ILD2 - p.F2 < 0 ))
                continue;
            end
            
            LB_E = p.E2 - F_F2 * ILC2;
            LB_E = max(0, LB_E);
            
            if any(LB_E(:) > p.x_max)
                continue;
            end
            
            F_E2 = LB_E + (p.x_max - LB_E) .* rand(p.dims.F_E2);
            
            pop(i).pos.F_F2 = F_F2; pop(i).vel.F_F2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_F2);
            pop(i).pos.F_E2 = F_E2; pop(i).vel.F_E2 = p.v_min + (p.v_max - p.v_min) * rand(p.dims.F_E2);
            is_EF2_valid = validate_EF(pop(i).pos, p, 2);
        end
        
        pop(i).pbest_is_valid = false; 
        pop(i) = evaluate_particle_fitness(pop(i), p, HUGE_PENALTY);
        pop(i) = update_individual_best(pop(i));
        [gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
            update_global_best(pop(i), gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty);
    end
end