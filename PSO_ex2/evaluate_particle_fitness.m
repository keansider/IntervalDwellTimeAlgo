function particle = evaluate_particle_fitness(particle, p, HUGE_PENALTY)  
    inputs = {particle.pos.F_B1, particle.pos.F_B2, ...
              particle.pos.F_F1, particle.pos.F_F2, ...
              particle.pos.F_A1, particle.pos.F_A2, ...
              particle.pos.F_E1, particle.pos.F_E2}; 
              
    [sol_main, err_main] = p.P_main(inputs);  
    if err_main == 0 && ~isnan(sol_main) && ~isinf(sol_main) 
        particle.is_valid = true;
        particle.fitness = sol_main;
        particle.penalty = 0;  
    else 
        particle.is_valid = false;
        particle.fitness = inf;   
        [sol_para, err_para] = p.P_para(inputs); 
        if err_para == 0 && ~isnan(sol_para) && ~isinf(sol_para)
            particle.penalty = sol_para; 
        else
            particle.penalty = HUGE_PENALTY;  
        end
    end 
end