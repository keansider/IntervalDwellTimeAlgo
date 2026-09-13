function result = find_opt_gamma(delta_xi, tmin, divide_number)
    mode_number = 2;
    P = build_yalmip_optimizer(tmin, divide_number, mode_number);
    
    xi_array = delta_xi : delta_xi : 1;
    gamma_array = NaN(size(xi_array));
    
    for k = 1 : length(xi_array)
        xi_val = xi_array(k);
        lambda_val = -log(xi_val) / tmin;
        
        [opt_g, errorcode] = P([lambda_val; xi_val]);
        
        if errorcode == 0
            gamma_array(k) = opt_g;
        end
    end
    
    result.xi = xi_array;
    result.gamma = gamma_array;
    
    [min_gamma, min_idx] = min(gamma_array);
    if ~isnan(min_gamma)
        min_xi = xi_array(min_idx);
        min_lambda = -log(min_xi) / tmin;
        fprintf('最优(最小) gamma = %.4f, 对应 xi = %.4f, lambda = %.4f\n', min_gamma, min_xi, min_lambda);
    else
        fprintf('未找到可行的 gamma 解\n');
    end
    
    figure;
    plot(xi_array, gamma_array, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    xlabel('\xi');
    ylabel('Optimal \gamma');
    grid on;
    
    current_ylim = ylim;
    if current_ylim(1) < 1
        ylim([current_ylim(1), 1]);
    end
end