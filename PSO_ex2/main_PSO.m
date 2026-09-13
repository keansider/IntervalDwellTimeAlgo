function main_PSO()
    clc
    clear all
    
    use_custom_first = 1;
    
    p = system_para(); 
    [p.P_main, p.P_para] = build_optimizers(p); 
    [pop, gbest_pos, gbest_fitness, gbest_penalty, gbest_is_valid] = init_population(p); 
    fields = fieldnames(p.dims);
    
    best_fitness_curve = zeros(1, p.max_iter);
    avg_fitness_curve = zeros(1, p.max_iter);
    best_penalty_curve = zeros(1, p.max_iter);
    avg_penalty_curve = zeros(1, p.max_iter);
    
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
        
        % 1. 记录最佳 Fitness (如果是 inf, MATLAB 画图会自动跳过)
        best_fitness_curve(iter) = gbest_fitness;
        
        % 2. 计算平均 Fitness (仅限有效个体)
        valid_idx = [pop.is_valid];
        if any(valid_idx)
            avg_fitness_curve(iter) = mean([pop(valid_idx).fitness]);
        else
            avg_fitness_curve(iter) = NaN; % 使用 NaN，画图时会自动留空
        end
        
        % 3. 计算最佳/平均 Penalty (仅限无效个体)
        invalid_idx = ~valid_idx;
        if any(invalid_idx)
            best_penalty_curve(iter) = min([pop(invalid_idx).penalty]);
            avg_penalty_curve(iter) = mean([pop(invalid_idx).penalty]);
        else
            best_penalty_curve(iter) = 0;
            avg_penalty_curve(iter) = 0;
        end
        
        print_iteration_stats(iter, gbest_is_valid, gbest_fitness, gbest_penalty, pop);
        
    end
    
    save('last_gbest_pos.mat', 'gbest_pos');
    
    timestamp_str = datestr(now, 'yyyymmdd_HHMMSS');
    data_filename = ['plot_data_', timestamp_str, '.mat'];
    save(data_filename, 'best_fitness_curve', 'avg_fitness_curve', 'best_penalty_curve', 'avg_penalty_curve');
    
    iters = 1:length(best_fitness_curve);

    %% 画图代码

    iters = 1:length(best_fitness_curve);

    fit_best_plot = best_fitness_curve;
    fit_best_plot(isinf(fit_best_plot) | fit_best_plot <= 0) = NaN;
    
    fit_avg_plot = avg_fitness_curve;
    fit_avg_plot(isinf(fit_avg_plot) | fit_avg_plot <= 0) = NaN;
    
    pen_best_plot = best_penalty_curve;
    pen_best_plot(isinf(pen_best_plot) | pen_best_plot <= 0) = NaN;
    
    pen_avg_plot = avg_penalty_curve;
    pen_avg_plot(isinf(pen_avg_plot) | pen_avg_plot <= 0) = NaN;
    
    invalid_fit = isnan(fit_best_plot);
    invalid_pen = isnan(pen_best_plot);
    
    figure;
    
    subplot(2, 1, 1);
    hold on;
    p1 = plot(iters, log(fit_best_plot), '-', 'Color', '#7A378A', 'LineWidth', 1);
    p2 = plot(iters, log(fit_avg_plot), '--', 'Color', '#E84446', 'LineWidth', 1);
    
    yl1 = [-1, 6];
    
    bg_y1 = yl1(1) * ones(size(iters));
    bg_y1(invalid_fit) = yl1(2);
    
    bg1 = bar(iters, bg_y1, 1, 'FaceColor', '#6A8EC9', 'EdgeColor', 'none', 'BaseValue', yl1(1), 'FaceAlpha', 0.3);
    uistack(bg1, 'bottom');
    
    ylim(yl1);
    xlim([0, 5000]);
    
    ylabel('ln(Fitness)');
    title('PSO Convergence Curves - Fitness');
    legend([p1, p2, bg1], {'Best Fitness', 'Average Fitness', 'Invalid Region'});
    grid on;
    box on;
    set(gca, 'Layer', 'top');
    
    subplot(2, 1, 2);
    hold on;
    p3 = plot(iters, log(pen_best_plot), '-', 'Color', '#7A378A', 'LineWidth', 1);
    p4 = plot(iters, log(pen_avg_plot), '--', 'Color', '#E84446', 'LineWidth', 1);
    
    yl2 = [-12, 5];
    
    bg_y2 = yl2(1) * ones(size(iters));
    bg_y2(invalid_pen) = yl2(2);
    
    bg2 = bar(iters, bg_y2, 1, 'FaceColor', '#6A8EC9', 'EdgeColor', 'none', 'BaseValue', yl2(1), 'FaceAlpha', 0.3);
    uistack(bg2, 'bottom');
    
    ylim(yl2);
    xlim([0, 5000]);
    
    xlabel('Iteration');
    ylabel('ln(Penalty)');
    title('PSO Convergence Curves - Penalty');
    legend([p3, p4, bg2], {'Best Penalty', 'Average Penalty', 'Invalid Region'});
    grid on;
    box on;
    set(gca, 'Layer', 'top');

    %% 输出代码
    
    fprintf('\n==================================================\n');
    if gbest_is_valid
        fprintf('           找到可行解！最优解矩阵如下：\n');
    else
        fprintf('  未找到可行解。以下为惩罚数最小(%g)的解矩阵：\n', gbest_penalty);
    end
    fprintf('==================================================\n'); 
    
    matrix_names = fieldnames(gbest_pos);
    for idx = 1:length(matrix_names)
        mat_name = matrix_names{idx};
        fprintf('\n矩阵 %s =\n', mat_name);
        disp(gbest_pos.(mat_name));
    end
    fprintf('==================================================\n\n');
end