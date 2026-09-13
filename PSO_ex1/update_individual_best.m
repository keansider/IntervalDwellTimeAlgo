function particle = update_individual_best(particle) 

    update_pbest = false; 
    if isempty(particle.pbest_is_valid)
        particle.pbest_is_valid = false;
    end
    if isempty(particle.pbest_fitness)
        particle.pbest_fitness = inf;
    end
    if isempty(particle.pbest_penalty)
        particle.pbest_penalty = inf;
    end
    
    if particle.is_valid  
        if ~particle.pbest_is_valid || (particle.fitness < particle.pbest_fitness)
            update_pbest = true;
        end
    else 
        if ~particle.pbest_is_valid && (particle.penalty < particle.pbest_penalty)
            update_pbest = true;
        end
    end 
    
    if update_pbest
        particle.pbest_pos = particle.pos;
        particle.pbest_is_valid = particle.is_valid;
        particle.pbest_fitness = particle.fitness;
        particle.pbest_penalty = particle.penalty;
    end
end