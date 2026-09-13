function is_valid = validate_EF(pos, p, index)
    if index == 1 
        aux1 = all(pos.F_E1(:) >= 0);
        aux2 = all(pos.F_F1(:) >= 0);
        m1 = pos.F_E1 - p.E1 + pos.F_F1 * (eye(2) - p.Lambda) * p.C1; 
        aux3 = all(m1(:) >= 0); 
        aux4 = all(abs(pos.F_E1(:)) >= 1e-4);
        aux5 = all(abs(pos.F_F1(:)) >= 1e-4);
        if aux1 == 1 && aux2 == 1 && aux3 == 1 && aux4 == 1 && aux5 == 1  
            is_valid = true; 
        else
            is_valid = false;
        end
    else
        aux1 = all(pos.F_E2(:) >= 0);
        aux2 = all(pos.F_F2(:) >= 0);
        m1 = pos.F_E2 - p.E2 + pos.F_F2 * (eye(2) - p.Lambda) * p.C2; 
        aux3 = all(m1(:) >= 0); 
        aux4 = all(abs(pos.F_E2(:)) >= 1e-4);
        aux5 = all(abs(pos.F_F2(:)) >= 1e-4);
        if aux1 == 1 && aux2 == 1 && aux3 == 1 && aux4 == 1 && aux5 == 1 
            is_valid = true; 
        else
            is_valid = false;
        end
    end
end