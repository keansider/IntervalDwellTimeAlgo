function print_iteration_stats(iter, gbest_is_valid, gbest_fitness, gbest_penalty, pop)
    valid_idx = [pop.is_valid];
    valid_fitnesses = [pop(valid_idx).fitness];
    invalid_penalties = [pop(~valid_idx).penalty];
    
    if isempty(valid_fitnesses)
        avg_fit_str = 'null';
    else
        avg_fit_str = sprintf('%.4f', mean(valid_fitnesses));
    end
    
    if isempty(invalid_penalties)
        avg_pen_str = 'null';
    else
        avg_pen_str = sprintf('%.4f', mean(invalid_penalties));
    end
    
    if gbest_is_valid
        best_fit_str = sprintf('%.4f', gbest_fitness);
        best_pen_str = 'null';
    else
        best_fit_str = 'null';
        best_pen_str = sprintf('%.4f', gbest_penalty);
    end
    
    fprintf('Iteration: %d, Best Fit=%s, Avg Fit=%s, Best Pen=%s, Avg Pen=%s\n', ...
        iter, best_fit_str, avg_fit_str, best_pen_str, avg_pen_str);
end