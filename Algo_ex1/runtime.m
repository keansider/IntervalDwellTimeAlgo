%% 算法运行时间大样本对比脚本 (简洁版)
% 清理工作区变量、命令窗口和历史绘图
clear; clc; close all;
%% 1. 参数设置
divide_number = 1;
tmin = 0.3;
delta_xi_list = [0.1, 0.01, 0.001];
num_trials = 100; % 对比次数
num_cases = length(delta_xi_list);
% 预分配内存存储运行时间
t_algo1 = zeros(num_cases, num_trials); 
%% 2. 开始大样本对比测试
fprintf('=== 开始算法耗时对比 (共 %d 组, 每组 %d 次) ===\n\n', num_cases, num_trials);
for i = 1:num_cases
    delta_xi = delta_xi_list(i);
    fprintf('>>> 开始测试参数: delta_xi = %g\n', delta_xi);
    
    for j = 1:num_trials
        % 运行算法 1
        res1 = algo_ex1(delta_xi, tmin, divide_number);
        t_algo1(i, j) = res1.t;
        % 完成算法1后立即打印
        % 修改后：
    fprintf('[%d/%d] delta_xi: %g | 用时: %.4f s\n', j, num_trials, delta_xi, res1.t);
    end
    
    fprintf('<<< 参数 delta_xi = %g 测试完毕！\n\n', delta_xi);
end
%% 3. 统计结果输出
fprintf('=== 实验运行完成，耗时统计结果 (单位: 秒) ===\n');
figure('Name', '算法运行时间对比', 'Position', [100, 100, 900, 400]);
for i = 1:num_cases
    delta_xi = delta_xi_list(i);
    
    % 计算平均时间
    mean_t1 = mean(t_algo1(i, :)); 
    
    % 打印最终结果
    fprintf('【delta_xi = %g】\n', delta_xi);
    fprintf('  算法1 (二分法)   平均耗时: %.4f s\n', mean_t1); 
    
    % 绘制箱线图
    subplot(1, num_cases, i);
    boxplot([t_algo1(i, :)']);
    title(sprintf('\\Delta\\xi = %g', delta_xi));
    ylabel('Run time (s)');
    grid on;
end
save(['plot_data_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat'], 't_algo1', 'delta_xi_list', 'num_cases', 'num_trials');