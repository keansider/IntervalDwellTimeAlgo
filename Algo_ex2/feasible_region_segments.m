%% 主代码

clear; 
clc; 
close all;

% 分段数

divide_numbers = [1, 2, 4, 9];


results_data = cell(length(divide_numbers), 1);

fprintf('开始计算不同分段下的可行域...\n');
fprintf('-----------------------------------------\n');

for i = 1:length(divide_numbers)
    div_num = divide_numbers(i);
    fprintf('正在计算分段数 divide_number = %d ...\n', div_num);
    
    % 初始化 

    tmin_list = [];
    tmax_list = [];
    
    % 设定初始参数

    tmin = 0.01;
    step = 0.002;
    delta_xi = 0.01; 
    
    while true

        % 调用你定义的 algo_ex1 函数

        res = algo_ex1(delta_xi, tmin, div_num);
        
        % 当无法找到满足 tmax > tmin 的可行解时，停止当前分段数的递增循环

        if res.flag == 0
            fprintf('  => [停止] tmin = %.3f 时未找到可行解 (tmax <= tmin)\n', tmin);
            break;
        end
        
        % 记录数据

        tmin_list(end+1) = tmin;
        tmax_list(end+1) = res.tmax;
        
        % 【新增】实时打印每一次计算的结果

        fprintf('  -> tmin = %.3f  |  计算得到 tmax = %.4f\n', tmin, res.tmax);
        
        % 递增 tmin

        tmin = tmin + step;
    end
    
    % 保存当前 divide_number 的结果

    results_data{i}.divide_number = div_num;
    results_data{i}.tmin_list = tmin_list;
    results_data{i}.tmax_list = tmax_list;
    
    fprintf('-----------------------------------------\n');
end

fprintf('计算完成！开始绘制可视化图表...\n');

%% 可视化部分 

figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on;
grid on;
box on;
colors = ['r', 'b', 'g', 'm'];
markers = ['o', 's', '^', 'd'];

for i = 1:length(divide_numbers)
    data = results_data{i};
    if ~isempty(data.tmin_list)
        plot(data.tmin_list, data.tmax_list, ...
             'Color', colors(i), ...
             'LineStyle', '-', ...
             'LineWidth', 1.5, ...
             'Marker', markers(i), ...
             'MarkerSize', 5, ...
             'MarkerFaceColor', colors(i), ...
             'DisplayName', sprintf('Segments = %d', data.divide_number));
    end
end

% 绘制一条辅助线：tmax = tmin （临界边界）

all_tmins = cellfun(@(x) x.tmin_list, results_data, 'UniformOutput', false);
if ~isempty(cell2mat(all_tmins))
    max_tmin = max(cell2mat(all_tmins));
    plot([0, max_tmin], [0, max_tmin], 'k--', 'LineWidth', 1.2, 'DisplayName', 't_{max} = t_{min}');
end

% 界面美化

xlabel('Minimum Bound t_{min}', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Maximum Time t_{max}', 'FontSize', 12, 'FontWeight', 'bold');
title('Feasible Region under Different Segmentations', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
set(gca, 'FontSize', 11);

hold off;