figure('Color', 'w', 'Position', [100, 100, 700, 500]);
hold on;

xlim([0, 1.8]);
ylim([0, 1.8]);

xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
for k = 1:length(xt)
    plot([xt(k) xt(k)], [0 1.8], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for k = 1:length(yt)
    plot([0 1.8], [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

colors = {[0.1, 0.1, 0.1], [0.4, 0.4, 0.4], [0.7, 0.7, 0.7], [0.9, 0.9, 0.9]};
num_slices = length(divide_numbers);
fill_handles = gobjects(1, num_slices);

for i = num_slices:-1:1
    data = results_data{i};
    tmin = data.tmin_list;
    tmax = data.tmax_list;
    
    if isempty(tmin)
        continue;
    end
    
    X_poly = [tmin, fliplr(tmin)];
    Y_poly = [tmax, fliplr(tmin)];
    
    % 绘制填充
    fill_handles(i) = fill(X_poly, Y_poly, colors{i}, 'FaceAlpha', 1, 'EdgeColor', 'none');
    
    % 绘制上下边界线
    plot(tmin, tmax, 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot(tmin, tmin, 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
end

xlabel('\delta_L', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('\delta_U', 'FontName', 'Times New Roman', 'FontSize', 12);

box on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Layer', 'top', 'LineWidth', 0.8);

valid_idx = arrayfun(@(x) isgraphics(x), fill_handles);
legend(fill_handles(valid_idx), strcat('N = ', string(divide_numbers(valid_idx))), ...
    'Location', 'northwest', 'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'on');

hold off;