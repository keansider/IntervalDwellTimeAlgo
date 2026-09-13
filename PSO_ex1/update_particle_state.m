function particle = update_particle_state(particle, gbest_pos, w, p) 
    
    MAX_RETRY = 20; 

    is_AB1_valid = false;
    retry = 0;
    temp_pos = particle.pos; 
    temp_vel = particle.vel;
    while ~is_AB1_valid && retry < MAX_RETRY
        [temp_vel, temp_pos] = update_single_matrix('F_A1', temp_vel, temp_pos, particle, gbest_pos, w, p);
        [temp_vel, temp_pos] = update_single_matrix('F_B1', temp_vel, temp_pos, particle, gbest_pos, w, p);
        is_AB1_valid = validate_AB(temp_pos, p, 1);
        retry = retry + 1;
    end 
    if is_AB1_valid
        particle.vel.F_A1 = temp_vel.F_A1; 
        particle.pos.F_A1 = temp_pos.F_A1;
        particle.vel.F_B1 = temp_vel.F_B1; 
        particle.pos.F_B1 = temp_pos.F_B1;
    else
        particle.vel.F_A1 = -rand(p.dims.F_A1) .* particle.vel.F_A1;
        particle.vel.F_B1 = -rand(p.dims.F_B1) .* particle.vel.F_B1;
        particle.pos.F_A1 = particle.pbest_pos.F_A1;
        particle.pos.F_B1 = particle.pbest_pos.F_B1;
    end

    is_AB2_valid = false;
    retry = 0;
    temp_pos = particle.pos; 
    temp_vel = particle.vel;
    while ~is_AB2_valid && retry < MAX_RETRY
        [temp_vel, temp_pos] = update_single_matrix('F_A2', temp_vel, temp_pos, particle, gbest_pos, w, p);
        [temp_vel, temp_pos] = update_single_matrix('F_B2', temp_vel, temp_pos, particle, gbest_pos, w, p);
        is_AB2_valid = validate_AB(temp_pos, p, 2);
        retry = retry + 1;
    end
    if is_AB2_valid
        particle.vel.F_A2 = temp_vel.F_A2; 
        particle.pos.F_A2 = temp_pos.F_A2;
        particle.vel.F_B2 = temp_vel.F_B2; 
        particle.pos.F_B2 = temp_pos.F_B2;
    else
        particle.vel.F_A2 = -rand(p.dims.F_A2) .* particle.vel.F_A2;
        particle.vel.F_B2 = -rand(p.dims.F_B2) .* particle.vel.F_B2;
        particle.pos.F_A2 = particle.pbest_pos.F_A2;
        particle.pos.F_B2 = particle.pbest_pos.F_B2;
    end

    is_EF1_valid = false;
    retry = 0;
    temp_pos = particle.pos; 
    temp_vel = particle.vel;
    while ~is_EF1_valid && retry < MAX_RETRY
        [temp_vel, temp_pos] = update_single_matrix('F_E1', temp_vel, temp_pos, particle, gbest_pos, w, p);
        [temp_vel, temp_pos] = update_single_matrix('F_F1', temp_vel, temp_pos, particle, gbest_pos, w, p);
        is_EF1_valid = validate_EF(temp_pos, p, 1);
        retry = retry + 1;
    end
    if is_EF1_valid
        particle.vel.F_E1 = temp_vel.F_E1; 
        particle.pos.F_E1 = temp_pos.F_E1;
        particle.vel.F_F1 = temp_vel.F_F1; 
        particle.pos.F_F1 = temp_pos.F_F1;
    else
        particle.vel.F_E1 = -rand(p.dims.F_E1) .* particle.vel.F_E1;
        particle.vel.F_F1 = -rand(p.dims.F_F1) .* particle.vel.F_F1;
        particle.pos.F_E1 = particle.pbest_pos.F_E1;
        particle.pos.F_F1 = particle.pbest_pos.F_F1;
    end
    
    is_EF2_valid = false;
    retry = 0;
    temp_pos = particle.pos; 
    temp_vel = particle.vel;
    while ~is_EF2_valid && retry < MAX_RETRY
        [temp_vel, temp_pos] = update_single_matrix('F_E2', temp_vel, temp_pos, particle, gbest_pos, w, p);
        [temp_vel, temp_pos] = update_single_matrix('F_F2', temp_vel, temp_pos, particle, gbest_pos, w, p);
        is_EF2_valid = validate_EF(temp_pos, p, 2);
        retry = retry + 1;
    end
    if is_EF2_valid
        particle.vel.F_E2 = temp_vel.F_E2; 
        particle.pos.F_E2 = temp_pos.F_E2;
        particle.vel.F_F2 = temp_vel.F_F2; 
        particle.pos.F_F2 = temp_pos.F_F2;
    else
        particle.vel.F_E2 = -rand(p.dims.F_E2) .* particle.vel.F_E2;
        particle.vel.F_F2 = -rand(p.dims.F_F2) .* particle.vel.F_F2;
        particle.pos.F_E2 = particle.pbest_pos.F_E2;
        particle.pos.F_F2 = particle.pbest_pos.F_F2;
    end
end