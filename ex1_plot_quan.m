%% 系统状态图

figure;

data_t = {out.sub1.time, out.sub2.time};
data_v = {squeeze(out.sub1.signals.values), squeeze(out.sub2.signals.values)};
sub_titles = {'subsystem1', 'subsystem2'};

for i = 1:2
    subplot(2, 2, i);
    hold on;
    
    h1 = plot(data_t{i}, data_v{i}(:, 1), '-', 'Color', '#59B78F', 'LineWidth', 1, 'DisplayName', '$x_1$');
    h2 = plot(data_t{i}, data_v{i}(:, 2), '-.', 'Color', '#6A8EC9', 'LineWidth', 1, 'DisplayName', '$x_2$');
    
    axis([0 15 0 3500]);
    
    xt = get(gca, 'XTick');
    yt = get(gca, 'YTick');
    xl = xlim;
    yl = ylim;
    
    for k = 1:length(xt)
        plot([xt(k) xt(k)], yl, '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    for k = 1:length(yt)
        plot(xl, [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    
    uistack([h1, h2], 'top');
    box on;
    set(gca, 'Layer', 'top', 'LineWidth', 0.8);
    title(sub_titles{i});
    xlabel('t/s');
    legend([h1, h2], 'Location', 'northwest', 'Interpreter', 'latex');
end

subplot(2, 2, [3, 4]);
hold on;

t = out.x.time;
v = squeeze(out.x.signals.values);
if size(v, 1) ~= length(t)
    v = v';
end

colors = {'#59B78F', '#6A8EC9', '#7A378A', '#E84446'};
styles = {'-', '--', ':', '-.'};
labels = {'$x_1$', '$x_2$', '$\tilde{x}_1$', '$\tilde{x}_2$'};
h = zeros(1, 4);

for j = 1:4
    h(j) = plot(t, v(:, j), 'Color', colors{j}, 'LineStyle', styles{j}, 'LineWidth', 1, 'DisplayName', labels{j});
end

axis([0 15 0 12]);

xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
xl = xlim;
yl = ylim;

for k = 1:length(xt)
    plot([xt(k) xt(k)], yl, '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for k = 1:length(yt)
    plot(xl, [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

uistack(h, 'top');
box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
xlabel('t/s');
legend('Location', 'northeast', 'Interpreter', 'latex');



%% 切换律


figure;
hold on;

t = out.law.time;
v = squeeze(out.law.signals.values);

h = plot(t, v, 'Color', '#E84446', 'LineWidth', 1, 'DisplayName', '$\sigma (t)$');

axis([0 15 0.5 3.5]);

xt = get(gca, 'XTick');
yt = get(gca, 'YTick');

for k = 1:length(xt)
    plot([xt(k) xt(k)], [0 3.5], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for k = 1:length(yt)
    plot([0 20], [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

uistack(h, 'top');

box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
title('Switching Law');
xlabel('t/s');
legend('Location', 'northeast', 'Interpreter', 'latex');




%% 输出差异

figure;
hold on;
t = out.error_z.time;
h1 = plot(t, out.error_z.signals.values(1,:), '-', 'Color', '#7A378A', 'LineWidth', 1, 'DisplayName', '$\tilde{z}_1(t)$');
h2 = plot(t, out.error_z.signals.values(2,:), '-.', 'Color', '#E84446', 'LineWidth', 1, 'DisplayName', '$\tilde{z}_2(t)$');
axis([0 15 0 140]);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
for k = 1:length(xt)
    plot([xt(k) xt(k)], [0 140], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for k = 1:length(yt)
    plot([0 15], [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
uistack([h1, h2], 'top');
box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
xlabel('t/s');
legend([h1, h2], 'Location', 'northeast', 'Interpreter', 'latex');





%% 量化输出

figure;

ty = out.y.time;
vy = squeeze(out.y.signals.values);
if size(vy, 1) ~= length(ty), vy = vy'; end

tqy = out.quan_y.time;
vqy = squeeze(out.quan_y.signals.values);
if size(vqy, 1) ~= length(tqy), vqy = vqy'; end

colors = {'#59B78F', '#6A8EC9', '#7A378A', '#E84446'};
styles = {'-', '-.'};

subplot(1, 2, 1);
hold on;
plot(ty, vy(:, 1), 'Color', colors{1}, 'LineStyle', styles{1}, 'LineWidth', 1, 'DisplayName', '$y_1(t)$');
plot(ty, vy(:, 2), 'Color', colors{2}, 'LineStyle', styles{2}, 'LineWidth', 1, 'DisplayName', '$y_2(t)$');
axis([0 15 0 13]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.3, 'Layer', 'top', 'LineWidth', 0.8);
box on;
xlabel('t/s');
ylabel('y');
title('Original Output Trajectories');
legend('Location', 'northeast', 'Interpreter', 'latex');

subplot(1, 2, 2);
hold on;
plot(tqy, vqy(:, 1), 'Color', colors{3}, 'LineStyle', styles{1}, 'LineWidth', 1, 'DisplayName', '$\tilde{y}_1(t)$');
plot(tqy, vqy(:, 2), 'Color', colors{4}, 'LineStyle', styles{2}, 'LineWidth', 1, 'DisplayName', '$\tilde{y}_2(t)$');
axis([0 15 0 13]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.3, 'Layer', 'top', 'LineWidth', 0.8);
box on;
xlabel('t/s');
ylabel('$\tilde{y}$', 'Interpreter', 'latex');
title('Quantized Output Trajectories');
legend('Location', 'northeast', 'Interpreter', 'latex');











