function main_PSO()
    clc
    clear all
    
    use_custom_first = 1;

    p = system_para();
    [p.P_main, p.P_para] = build_optimizers(p);
    
    [pop, gbest_pos, gbest_fitness, gbest_penalty, gbest_is_valid] = init_population(p); 
    fields = fieldnames(p.dims);
    num_fields = length(fields);
    convergence_curve = zeros(1, p.max_iter);
    
    HUGE_PENALTY = 1e8;
    
    [pop, gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
        inject_custom_particle(pop, use_custom_first, p, HUGE_PENALTY, gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty);

    for iter = 1 : p.max_iter
        w = p.w_max - (p.w_max - p.w_min) * (iter / p.max_iter);
        
        for i = 1:p.pop_size
            pop(i) = update_particle_state(pop(i), gbest_pos, w, p); 
            pop(i) = evaluate_particle_fitness(pop(i), p, HUGE_PENALTY);
            pop(i) = update_individual_best(pop(i));
            [gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
                update_global_best(pop(i), gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty);
        end
        
        convergence_curve(iter) = gbest_fitness;
        print_iteration_stats(iter, gbest_is_valid, gbest_fitness, gbest_penalty, pop);
    end
    
    save('last_gbest_pos.mat', 'gbest_pos');

    figure;
    plot(convergence_curve, 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Best Fitness');
    title('PSO Convergence Curve');
    grid on;

    if gbest_is_valid
        fprintf('\n==================================================\n');
        fprintf('           找到可行解！最优解矩阵如下：\n');
        fprintf('==================================================\n'); 
        
        matrix_names = fieldnames(gbest_pos);
        for idx = 1:length(matrix_names)
            mat_name = matrix_names{idx};
            fprintf('\n矩阵 %s =\n', mat_name);
            disp(gbest_pos.(mat_name));
        end
        fprintf('==================================================\n\n');
    else
        fprintf('\n[提示] 迭代结束，未找到满足全部约束的可行解。\n\n');
    end
end