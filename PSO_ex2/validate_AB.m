function is_valid = validate_AB(pos, p, index)
    if index == 1 
        aux1 = double(all(pos.F_A1(~eye(size(pos.F_A1))) >= 0));
        aux2 = all(pos.F_B1(:) >= 0);
        m1 = pos.F_A1 - p.A1 + pos.F_B1 * (eye(3) - p.Lambda) * p.C1;
        m2 = pos.F_B1 * (eye(3) - p.Lambda) * p.D1 - p.B1;
        aux3 = all(m1(:) >= 0);
        aux4 = all(m2(:) >= 0);
        aux5 = all(abs(pos.F_A1(:)) >= 1e-4);
        aux6 = all(abs(pos.F_B1(:)) >= 1e-4);
        if aux1 == 1 && aux2 == 1 && aux3 == 1 && aux4 == 1 && aux5 == 1 && aux6 == 1
            is_valid = true; 
        else
            is_valid = false;
        end
    else
        aux1 = double(all(pos.F_A2(~eye(size(pos.F_A2))) >= 0));
        aux2 = all(pos.F_B2(:) >= 0);
        m1 = pos.F_A2 - p.A2 + pos.F_B2 * (eye(3) - p.Lambda) * p.C2;
        m2 = pos.F_B2 * (eye(3) - p.Lambda) * p.D2 - p.B2;
        aux3 = all(m1(:) >= 0);
        aux4 = all(m2(:) >= 0);
        aux5 = all(abs(pos.F_A2(:)) >= 1e-4);
        aux6 = all(abs(pos.F_B2(:)) >= 1e-4);
        if aux1 == 1 && aux2 == 1 && aux3 == 1 && aux4 == 1 && aux5 == 1 && aux6 == 1
            is_valid = true; 
        else
            is_valid = false;
        end
    end
end