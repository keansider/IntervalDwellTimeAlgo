function [gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty] = ...
    update_global_best(particle, gbest_pos, gbest_is_valid, gbest_fitness, gbest_penalty)

    update_gbest = false;
    if particle.is_valid 
        if ~gbest_is_valid || (particle.fitness < gbest_fitness)
            update_gbest = true;
        end
    else 
        if ~gbest_is_valid && (particle.penalty < gbest_penalty)
            update_gbest = true;
        end
    end 
    if update_gbest
        gbest_pos = particle.pos;
        gbest_is_valid = particle.is_valid;
        gbest_fitness = particle.fitness;
        gbest_penalty = particle.penalty;
    end
end