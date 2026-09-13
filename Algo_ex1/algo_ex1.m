function result = algo_ex1(delta_xi, tmin, divide_number)

    %% 时间记录

    t1 = tic;

    %% 参数初始化

    delta_lambda = 0.0001;
    epslion = 0.0001;
    mode_number = 2;
    xi = 1;
    lambda_left = 0;
    lambda_right = (floor(-(log(delta_xi) / (tmin * delta_lambda))) + 1) * delta_lambda;
    N = round(1 / delta_xi); 
    best_record = [-inf; -inf; -inf]; 
    
    P = build_yalmip_optimizer(tmin, divide_number, mode_number);
    
    %% 辅助矩阵初始化

    F = zeros(1, N - 1);
    E = -inf(3, N - 1); 
    
    for b = 1 : N - 1
        xi = xi - delta_xi; 
        l_left = lambda_left;
        if b ~= 1 && F(b - 1) > 0 
            l_left = F(b - 1); 
        end
        l_right = lambda_right; 
        best_valid_lambda = -1;  
        
        %% 二分法搜索

        while (l_right - l_left) > epslion
            mid_lambda = 0.5 * (l_left + l_right);
            [~, errorcode] = P([mid_lambda; xi]);
            if errorcode == 0
                best_valid_lambda = mid_lambda;  
                l_right = mid_lambda;       
            else
                l_left = mid_lambda;  
            end
        end 
        
        %% 记录数据

        if best_valid_lambda ~= -1 
            F(b) = best_valid_lambda;
            E(1, b) = xi;
            E(2, b) = best_valid_lambda;
            E(3, b) = -log(xi) / best_valid_lambda;
        else 
            if b > 1
                F(b) = F(b - 1);
            else
                F(b) = 0; 
            end 
        end 
        
        if E(3, b) >= best_record(3)
            best_record = E(1:3, b);
        end
    end

    runtime = toc(t1);
    result.t = runtime;
    
    %% 结果判定

    if - log(best_record(1)) / best_record(2) <= tmin
        result.flag = 0; 
    else
        result.xi = best_record(1);
        result.lambda = best_record(2);
        result.tmax = - log(best_record(1)) ./ best_record(2);
        result.flag = 1;
        result.E = E;
    end
end