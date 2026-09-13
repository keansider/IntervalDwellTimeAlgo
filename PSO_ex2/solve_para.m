function [result, flag] = solve_para(lambda, tmin, xi, divide_number, mode_number, filter_matrix)

% clc
% clear   


%% 时间记录

% t1 = tic;

%% 辅助矩阵

w = ones(6, 1);
r = ones(3, 1);
f_vector = [];
I_y = eye(3);
O = zeros(3);


%% 系统参数

A{1} = [-3.60 1.02 0.52; 0.97 -3.40 0.82; 0.47 0.78 -3.20];
A{2} = [-3.10 1.02 0; 0.97 -3.40 0.82; 0 0.78 -2.70];

B{1} = [0.1032 0.2269 0.2159]';
B{2} = [0.1681 0.1895 0.1376]';

C{1} = [0.7165 0.8395 0.7548; 0.7266 0.7146 0.8842; 0.8848 0.7452 0.7991];
C{2} = [0.8074 0.8001 0.7065; 0.8424 0.7494 0.7586; 0.7303 0.8909 0.8349];

D{1} = [0.1798 0.1507 1.1006]';
D{2} = [0.1442 0.1069 0.2358]';

E{1} = [0.6351 0.4383 0.5213; 0.6515 0.6732 0.6919; 0.4194 0.6849 0.5771];
E{2} = [0.6205 0.4271 0.6223; 0.5871 0.6387 0.6665; 0.6948 0.4019 0.6692];

F{1} = [0.1835 0.2158 0.1677]';
F{2} = [0.1457 0.1087 0.2149]';

Lambda = diag([0.26 0.15 0.17]);


%% 测试参数

% xi = 0.72;
% lambda = 0.9133;
% tmin = 0.7;
% divide_number = 1;
% mode_number = 2;  


%% 滤波器参数

F_B{1} = filter_matrix.F_B1;
F_B{2} = filter_matrix.F_B2;

F_F{1} = filter_matrix.F_F1;
F_F{2} = filter_matrix.F_F2;

F_A{1} = filter_matrix.F_A1;
F_A{2} = filter_matrix.F_A2;

F_E{1} = filter_matrix.F_E1;
F_E{2} = filter_matrix.F_E2;

% F_B{1} = [5.2137 0.7149 2.2006; 3.4607 1.7133 3.0770; 5.2137 5.1955 4.7216];
% F_B{2} = [1.7622 3.6445 4.5118; 2.5524 2.7137 5.2137; 0.3997 3.4269 5.2137];
% F_A{1} = [-5.9982 0.8052 0.7690; 0.4206 -5.7469 0.0011; 0.0945 0.6742 -6.9469];
% F_A{2} = [-6.2146 1.6123 1.6541; 0.0323 -8.4346 0.1973; 0.2436 0.2137 -8.1198];
% F_F{1} = [4.3313 0.5880 0.0028; 5.2137 3.6339 2.5655; 2.0859 0.5744 0.5153];
% F_F{2} = [2.6931 4.0988 2.8421; 2.7260 2.7137 5.2137; 5.1253 0.8383 0.3740];
% F_E{1} = [0.8655 2.8305 0.8383; 2.6806 3.6871 0.1271; 0.0351 1.2407 3.5741];
% F_E{2} = [2.5967 1.8164 1.4696; 1.2193 0.3961 0.0848; 0.0067 1.0066 2.5728];

%% 矩阵参数设置

hat_A = cell(mode_number); 
hat_B = cell(mode_number); 
hat_C = cell(mode_number); 
hat_D = cell(mode_number); 

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
        f{i, j} = sdpvar(6, 1);
    end
end

gamma = sdpvar(1);
lar = sdpvar(1);

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
        Lvi = hat_A{i}' * f{i, q} + aux1 - lambda * f{i, q} + hat_C{i}' * r - lar * w;
        Cons{i, q} = Lvi <= -e1;
    end
end

Constraints = [Constraints, Cons{:}];

%% (2).2式 # 5-6

Cons = cell(mode_number, divide_number);

for i = 1 : mode_number
    for q = 1 : divide_number
        aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
        Lvi = hat_A{i}' * f{i, q + 1} + aux1 - lambda * f{i, q + 1} + hat_C{i}' * r - lar * w;
        Cons{i, q} = Lvi <= -e1; 
    end
end

Constraints = [Constraints, Cons{:}];

%% (2)式 # 7-8

Cons = cell(1, mode_number);

for i = 1 : mode_number 
    Lvi = hat_A{i}' * f{i, divide_number + 1} - lambda * f{i, divide_number + 1} + hat_C{i}' * r - lar * w;
    Cons{1, i} = Lvi <= -e1; 
end

Constraints = [Constraints, Cons{:}];

%% (3)式 # 9-

Cons = cell(mode_number, divide_number + 1);

for i = 1 : mode_number 
    for q = 1 : divide_number + 1
        Lvi = hat_B{i}' * f{i, q} - gamma * 1 + hat_D{i}' * r - lar * r;
        Cons{i, q} = Lvi <= -e1; 
    end
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

Constraints = [Constraints, gamma >= e1];

%% 目标函数
 
Objective = lar;


%% 求解器设置 

options = sdpsettings('solver', 'mosek', ...
                      'cachesolvers', 1, ...        % 【新增】让 YALMIP 缓存已找到的 MOSEK 求解器，避免每次迭代重新检索
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
    result.aux = value(lar);
    for i = 1: mode_number
        for j = 1: divide_number + 1
            f_vector = [f_vector; value(f{i, j})];
        end
    end
    result.f_vector = f_vector;
elseif sol.problem ~= 0
    result.info = ['求解失败，YALMIP 报错编号：' num2str(sol.problem) ', 原因：' sol.info];
    flag = 0;
end
 
% runtime = toc(t1);
% fprintf('代码运行时间：%.4f 秒\n', runtime);



end







