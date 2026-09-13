function [pop, gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = inject_custom_particle(pop, flag, p, HUGE_PENALTY, gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty)
    if flag == 1 && exist('last_gbest_pos.mat', 'file')
        data = load('last_gbest_pos.mat');
        pop(1).pos = data.gbest_pos;
        pop(1).pbest_is_valid = false;
        pop(1) = evaluate_particle_fitness(pop(1), p, HUGE_PENALTY);
        pop(1) = update_individual_best(pop(1));
        [gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
            update_global_best(pop(1), gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty);
    end
end