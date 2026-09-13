%% 3D 可视化部分 (保留 XYZ 主轴线，去除外围框线)

% 创建画布并设置白色背景
figure('Color', 'w', 'Position', [100, 100, 900, 700]);
hold on; 
grid off; % 关闭网格线，使画面更干净
view(45, 30); % 设置初始的 3D 观察视角

num_slices = length(divide_numbers);

% 设置颜色映射
cmap = colormap('parula'); 
colors = cmap(round(linspace(1, 256, num_slices)), :);

% 直接使用真实的 divide_numbers 作为 Y 轴真实坐标
y_positions = divide_numbers;

% --- 1. 绘制带有上下界的多边形填充面 ---
for i = 1:num_slices
    data = results_data{i};
    tmin = data.tmin_list;
    tmax = data.tmax_list;
    
    if isempty(tmin)
        continue;
    end
    
    y_val = y_positions(i);
    
    % 构建闭合多边形坐标 (上界为 tmax，下界为 tmin)
    X_poly = [tmin, fliplr(tmin)];
    Y_poly = repmat(y_val, 1, 2*length(tmin));
    Z_poly = [tmax, fliplr(tmin)];
    
    % 绘制半透明的切片填充面
    fill3(X_poly, Y_poly, Z_poly, colors(i,:), ...
        'FaceAlpha', 0.45, 'EdgeColor', 'none');
    
    % 绘制顶部边界实线 (tmax) - 粗细设为 1
    plot3(tmin, repmat(y_val, 1, length(tmin)), tmax, ...
        'Color', colors(i,:), 'LineWidth', 1.5);
    
    % 绘制底部边界实线 (tmin) - 粗细设为 1
    plot3(tmin, repmat(y_val, 1, length(tmin)), tmin, ...
        'Color', colors(i,:), 'LineWidth', 1.5, 'LineStyle', '-.');
end

% --- 2. 预处理数据，跨块绘制连线 ---
% 找到所有数据中 tmin 的最大值
max_tmin_global = 0;
for i = 1:num_slices
    if ~isempty(results_data{i}.tmin_list)
        max_tmin_global = max(max_tmin_global, max(results_data{i}.tmin_list));
    end
end

% 步长为 0.3
target_x_vals = 0.3:0.3:max_tmin_global;

for j = 1:length(target_x_vals)
    x_target = target_x_vals(j);
    
    cross_x = [];
    cross_y = [];
    cross_z = [];
    
    % 在每个块中寻找是否存在当前的 x_target
    for i = 1:num_slices
        tmin = results_data{i}.tmin_list;
        tmax = results_data{i}.tmax_list;
        
        % 浮点数容差匹配
        idx = find(abs(tmin - x_target) < 1e-4);
        if ~isempty(idx)
            cross_x(end+1) = tmin(idx(1));
            cross_y(end+1) = y_positions(i);
            cross_z(end+1) = tmax(idx(1));
        end
    end
    
    % 只有当包含该数据点的分段数大于 1 时，才绘制折线和标签
    if length(cross_x) > 1
        % 绘制贯穿不同块的连线 - 粗细设为 1，标记颜色换为 #20c1b1
        plot3(cross_x, cross_y, cross_z, '--o', ...
            'Color', [0.3, 0.3, 0.3], ... 
            'LineWidth', 1.5, ...
            'MarkerSize', 5, ...
            'MarkerEdgeColor', '#20c1b1', ...    
            'MarkerFaceColor', '#20c1b1');       
        
        % 在节点上添加数据标签 - 字体颜色换为黑色 'k'
        for k = 1:length(cross_x)
            text(cross_x(k), cross_y(k), cross_z(k) + 0.05, sprintf('%.4f', cross_z(k)), ...
                'Color', 'k', 'HorizontalAlignment', 'center');
        end
    end
end

% --- 3. 界面美化与坐标轴设置 ---
xlabel('\delta_L');
ylabel('Divide Number');
zlabel('\delta_U');
title('Maximum Dwell Time Bounds under Different Segments');

% 强制设置坐标轴的范围
xlim([0, 1.8]);
ylim([1, 9]);  
zlim([0, 1.8]);

% 为了确保 Y 轴刻度只在你的分段数上显示
yticks(divide_numbers);
yticklabels(string(divide_numbers));

% 添加并配置右侧的 Colorbar
cb = colorbar;
caxis([min(divide_numbers), max(divide_numbers)]);
set(cb, 'Ticks', divide_numbers, 'TickLabels', string(divide_numbers));

% 设置图表边框和透视视角
set(gca, 'Box', 'off', ... % <--- 【核心修改】：关闭包围盒，只保留 XYZ 轴的线
         'Color', [0.96, 0.96, 0.96], ... % 坐标轴墙壁背景保留浅灰色
         'Projection', 'perspective');
hold off;