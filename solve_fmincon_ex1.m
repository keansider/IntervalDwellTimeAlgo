function [result, flag] = solve_fmincon_ex1(lambda, tmin, xi, divide_number, mode_number, filter_matrix)

% clc
% clear   

%% 时间记录

t1 = tic;

%% 辅助矩阵

w = ones(2, 1);
r = ones(2, 1);
f_vector = [];
I_y = eye(2);
O = zeros(2);


%% 系统参数

A{1} = [-2.487 0.231; 0.322 0.402];
A{2} = [0.366 0.311; 0.392 -2.137];

B{1} = [0.143 0.223]';
B{2} = [0.147 0.303]';

C{1} = [0.758 0.826; 0.822 0.874];
C{2} = [0.898 0.842; 0.762 0.827];

D{1} = [0.145 0.269]';
D{2} = [0.125 0.249]';

E{1} = [0.658 0.682; 0.694 0.427];
E{2} = [0.648 0.462; 0.623 0.617];

F{1} = [0.225 0.261]';
F{2} = [0.185 0.279]';

Lambda = diag([0.1, 0.2]);


%% 测试参数

xi = 0.7259;
lambda = 0.9123;
tmin = 0.3;
divide_number = 1;
mode_number = 2;  


%% 滤波器参数

% F_B{1} = filter_matrix.F_B1;
% F_B{2} = filter_matrix.F_B2;
% 
% F_F{1} = filter_matrix.F_F1;
% F_F{2} = filter_matrix.F_F2;
% 
% F_A{1} = filter_matrix.F_A1;
% F_A{2} = filter_matrix.F_A2;
% 
% F_E{1} = filter_matrix.F_E1;
% F_E{2} = filter_matrix.F_E2; 

F_B{1} = [3.0160 0.0029; 0.5439 0.0256];
F_B{2} = [1.4346 0.0303; 2.7971 0.0121];
F_A{1} = [-4.5377 0.0035; 0.0270 -0.0168];
F_A{2} = [-0.8108 0.0059; 0.0103 -4.2193];
F_F{1} = [3.8209 3.8209; 3.7631 3.8209];
F_F{2} = [1.3228 3.0085; 3.8209 3.8209];
F_E{1} = [2.5763 3.8209; 3.8209 1.1209];
F_E{2} = [3.8209 3.8209; 0.8349 3.8209];

%% 矩阵参数设置

for i = 1 : mode_number
    hat_A{i} = [A{i} O; F_A{i} - A{i} + F_B{i} * (I_y + Lambda) * C{i} F_A{i}];
    hat_B{i} = [B{i}; F_B{i} * (I_y + Lambda) * D{i} - B{i}];
    hat_C{i} = [F_E{i} - E{i} + F_F{i} * (I_y + Lambda) * C{i} F_E{i}];
    hat_D{i} = [F_F{i} * (I_y + Lambda) * D{i} - F{i}];
end


%% 预设量

Constraints = [];


%% 精度参数

e1 = 1e-4;
e2 = 1e-8;


%% 创建变量

f = cell(mode_number, divide_number + 1);  

for i = 1 : mode_number
    for j = 1 : divide_number + 1
        f{i, j} = sdpvar(4, 1);
    end
end


%% (1)式 # 1-2

Cons = cell(mode_number, mode_number);

for i = 1 : mode_number
    for j = 1 : mode_number
        if i ~= j 
            Lvi = f{j, 1} - xi * f{i, divide_number + 1};
            Cons{i, j} = Lvi <= -e1;
        end
    end
end

Constraints = [Constraints, Cons{:}];

%% (2).1式 # 3-4

Cons = cell(mode_number, divide_number);

for i = 1 : mode_number
    for q = 1 : divide_number
        aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
        Lvi = hat_A{i}' * f{i, q} + aux1 - lambda * f{i, q};
        Cons{i, q} = Lvi <= -e1;
    end
end

Constraints = [Constraints, Cons{:}];

%% (2).2式 # 5-6

Cons = cell(mode_number, divide_number);

for i = 1 : mode_number
    for q = 1 : divide_number
        aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
        Lvi = hat_A{i}' * f{i, q + 1} + aux1 - lambda * f{i, q + 1};
        Cons{i, q} = Lvi <= -e1; 
    end
end

Constraints = [Constraints, Cons{:}];

%% (2)式 # 7-8

Cons = cell(1, mode_number);

for i = 1 : mode_number 
    Lvi = hat_A{i}' * f{i, divide_number + 1} - lambda * f{i, divide_number + 1};
    Cons{1, i} = Lvi <= -e1; 
end

Constraints = [Constraints, Cons{:}];

%% 正约束

Cons = cell(mode_number, divide_number + 1);

for i = 1 : mode_number
    for j = 1 : divide_number + 1 
        Lvi = f{i, j};
        Cons{i, j} = Lvi >= e1; 
    end
end

Constraints = [Constraints, Cons{:}];
 
%% 辅助约束

Cons = cell(mode_number, divide_number + 1);

for i = 1 : mode_number
    for j = 1 : divide_number + 1 
        Lvi = f{i, j};
        Cons{i, j} = Lvi <= 1000; 
    end
end

Constraints = [Constraints, Cons{:}];


%% 目标函数
 
Objective = 0;


%% 求解器设置

options = sdpsettings('solver', 'mosek', ...
                      'verbose', 0, ...             % 关闭 YALMIP 表层打印
                      'debug', 0, ...               % 关闭 YALMIP 调试信息
                      'showprogress', 0, ...        % 关闭进度条提示
                      'mosek.MSK_IPAR_LOG', 0, ...  % 【关键】彻底关闭 MOSEK 的底层迭代日志打印
                      'mosek.MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-7,...
                      'mosek.MSK_DPAR_INTPNT_TOL_PFEAS', 1e-7,...
                      'mosek.MSK_DPAR_INTPNT_TOL_DFEAS', 1e-7);


%% 求解

sol = optimize(Constraints, Objective, options);
% check(Constraints)


%% 结果输出

if sol.problem == 0  
    result.success = true;
    result.info = '求解成功';
    flag = 1;
    for i = 1: mode_number
        for j = 1: divide_number + 1
            f_vector = [f_vector; value(f{i, j})];
        end
    end
    result.f_vector = f_vector;
    f_vector
elseif sol.problem ~= 0
    result.info = ['求解失败，YALMIP 报错编号：' num2str(sol.problem) ', 原因：' sol.info];
    flag = 0;
elseif any(res < -1e-4)
    result.info = '求解器返回成功，但部分约束未严格满足';
    flag = 0;
end

if sol.problem == 0
    disp('✅ 可行解成功找到！')
elseif sol.problem == 1
    disp('❌ 求解器无法判断可行性（可能发散）')
else
    disp(['⚠️ 求解失败，原因编号：', num2str(sol.problem)])
end

runtime = toc(t1);
fprintf('代码运行时间：%.4f 秒\n', runtime);



end







