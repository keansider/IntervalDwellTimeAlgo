function [temp_vel, temp_pos] = update_single_matrix(f_name, temp_vel, temp_pos, particle, gbest_pos, w, p)
    current_dim = p.dims.(f_name); 
    
    r1 = rand(current_dim); 
    r2 = rand(current_dim);
    
    % 更新速度

    aux1 = w * particle.vel.(f_name);
    aux2 = p.c1 * r1 .* (particle.pbest_pos.(f_name) - particle.pos.(f_name));
    aux3 = p.c2 * r2 .* (gbest_pos.(f_name) - particle.pos.(f_name));
    temp_vel.(f_name) = aux1 + aux2 + aux3;
        
    % 防止越界
    
    temp_vel.(f_name) = max(min(temp_vel.(f_name), p.v_max), p.v_min);
    temp_pos.(f_name) = particle.pos.(f_name) + temp_vel.(f_name);
    temp_pos.(f_name) = max(min(temp_pos.(f_name), p.x_max), p.x_min);
    
end