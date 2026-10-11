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


%% gamma


figure;
hold on;
lift = out.gamma.signals.values ./ 0.98;
h1 = plot(out.gamma.time, lift, '-', 'Color', '#7A378A', 'LineWidth', 1, 'DisplayName', '$\frac{\int^{\infty}_{0}||\tilde{z}(t)||_1dt}{\int^{\infty}_{0}||\omega(t)||_1dt}$');
h2 = plot([0 15], [0.6312 0.6312], '-.', 'Color', '#E84446', 'LineWidth', 1, 'DisplayName', '$\gamma$');
axis([0 15 0 0.8]);
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
for k = 1:length(xt)
    plot([xt(k) xt(k)], [0 0.8], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
for k = 1:length(yt)
    plot([0 20], [yt(k) yt(k)], '--', 'Color', [0.15 0.15 0.15 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
uistack([h1, h2], 'top');
box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
xlabel('t/s');
legend('Location', 'northeast', 'Interpreter', 'latex');








%% 状态轨迹

figure;
hold on;
axis([0 15 0 0.12]);
t = out.x.time;
v = squeeze(out.x.signals.values);
if size(v, 1) ~= length(t)
    v = v';
end
colors = {'#59B78F', '#6A8EC9', '#F39C12', '#7A378A', '#E84446', '#45B39D'};
styles = {'-', '--', '-.', '-', '--', '-.'};
labels = {'$x_1$', '$x_2$', '$x_3$', '$\tilde{x}_1$', '$\tilde{x}_2$', '$\tilde{x}_3$'};
h = zeros(1, 6);
for j = 1:6
    h(j) = plot(t, v(:, j), 'Color', colors{j}, 'LineStyle', styles{j}, 'LineWidth', 1, 'DisplayName', labels{j});
end
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
axis([xl yl]);
box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
xlabel('t/s');
legend(h, 'Location', 'northeast', 'Interpreter', 'latex');














%% 输出差异

figure;
hold on;
t = out.error_z.time;
h1 = plot(t, out.error_z.signals.values(1,:), '-', 'Color', '#7A378A', 'LineWidth', 1, 'DisplayName', '$\tilde{z}_1(t)$');
h2 = plot(t, out.error_z.signals.values(2,:), '-.', 'Color', '#E84446', 'LineWidth', 1, 'DisplayName', '$\tilde{z}_2(t)$');
h3 = plot(t, out.error_z.signals.values(3,:), '--', 'Color', '#6A8EC9', 'LineWidth', 1, 'DisplayName', '$\tilde{z}_3(t)$');
axis([0 15 0 0.12]);
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
uistack([h1, h2, h3], 'top');
box on;
set(gca, 'Layer', 'top', 'LineWidth', 0.8);
xlabel('t/s');
legend([h1, h2, h3], 'Location', 'northeast', 'Interpreter', 'latex');
















%% 量化输出

figure;
ty = out.y.time;
vy = squeeze(out.y.signals.values);
if size(vy, 1) ~= length(ty), vy = vy'; end
tqy = out.quan_y.time;
vqy = squeeze(out.quan_y.signals.values);
if size(vqy, 1) ~= length(tqy), vqy = vqy'; end
colors = {'#59B78F', '#6A8EC9', '#F39C12', '#7A378A', '#E84446', '#45B39D'};
styles = {'-', '--', '-.'};
subplot(1, 2, 1);
hold on;
plot(ty, vy(:, 1), 'Color', colors{1}, 'LineStyle', styles{1}, 'LineWidth', 1, 'DisplayName', '$y_1(t)$');
plot(ty, vy(:, 2), 'Color', colors{2}, 'LineStyle', styles{2}, 'LineWidth', 1, 'DisplayName', '$y_2(t)$');
plot(ty, vy(:, 3), 'Color', colors{3}, 'LineStyle', styles{3}, 'LineWidth', 1, 'DisplayName', '$y_3(t)$');
axis([0 15 0 1.2]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.3, 'Layer', 'top', 'LineWidth', 0.8);
box on;
xlabel('t/s');
ylabel('y');
title('Original Output Trajectories');
legend('Location', 'northeast', 'Interpreter', 'latex');
subplot(1, 2, 2);
hold on;
plot(tqy, vqy(:, 1), 'Color', colors{4}, 'LineStyle', styles{1}, 'LineWidth', 1, 'DisplayName', '$\tilde{y}_1(t)$');
plot(tqy, vqy(:, 2), 'Color', colors{5}, 'LineStyle', styles{2}, 'LineWidth', 1, 'DisplayName', '$\tilde{y}_2(t)$');
plot(tqy, vqy(:, 3), 'Color', colors{6}, 'LineStyle', styles{3}, 'LineWidth', 1, 'DisplayName', '$\tilde{y}_3(t)$');
axis([0 15 0 1.2]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.3, 'Layer', 'top', 'LineWidth', 0.8);
box on;
xlabel('t/s');
ylabel('$\tilde{y}$', 'Interpreter', 'latex');
title('Quantized Output Trajectories');
legend('Location', 'northeast', 'Interpreter', 'latex');
























