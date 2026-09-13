function result = algo_ex2(delta_xi, tmin, divide_number)

    t1 = tic;

    %% 参数初始化

    delta_lambda = 0.001;
    epslion = 0.0001;
    mode_number = 2;
    xi = 1;
    lambda_left = 0;
    lambda_right = (round(-(log(delta_xi) / (tmin * delta_lambda))) + 1) * delta_lambda;
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
        
        %% 二分法搜索 (此处已根据要求微调为：暴力迭代搜索)
        
        current_lambda = delta_lambda; % 从 0.0001 开始
        
        while current_lambda <= l_right
            [~, errorcode] = P([current_lambda; xi]);
            if errorcode == 0
                best_valid_lambda = current_lambda;  
                break; % 因为是从小往大找，找到的第一个有效值即为"最小允许值"，直接跳出循环 
            else
                current_lambda = current_lambda + delta_lambda; % 步长 0.0001 往上递增
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
    end
end