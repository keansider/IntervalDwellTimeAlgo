function [result, flag] = solve_fmincon_ex2(lambda, tmin, xi, divide_number, mode_number, filter_matrix)

% clc
% clear   

%% 时间记录

t1 = tic;

%% 辅助矩阵

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

xi = 0.98;
lambda = 0.0144;
tmin = 1.4;
divide_number = 100;
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

F_B{1} = [0.0786 1.6002 0.0003; 0.1723 1.3878 0.0287; 0.0067 1.3721 0.0430];
F_A{1} = [-4.6295 0.0071 0.0008; 0.0016 -4.3654 0.0044; 0.0001 0.0001 -4.2634];
F_B{2} = [0.0175 0.2783 3.2173; 0.2456 2.6515 0.0064; 0.0003 0.1476 0.9268];
F_A{2} = [-5.2349 0.0008 0.0019; 0.0033 -5.2349 0.0011; 0.0001 0.0005 -3.4376];
F_F{1} = [0.0001 0.9472 0.0680; 0.0001 1.0120 0.0943; 0.2159 0.8743 0.0295];
F_E{1} = [0.0001 0.0001 0.0001; 0.0001 0.0001 0.0001; 0.0001 0.0014 0.0001];
F_F{2} = [0.0001 0.3889 0.5638; 0.0001 0.6742 0.2827; 0.0001 0.0667 1.0670];
F_E{2} = [0.0002 0.0001 0.0001; 0.0001 0.0001 0.0358; 0.0002 0.0001 0.0001];

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
        f{i, j} = sdpvar(6, 1);
    end
end

gamma = sdpvar(1);

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
        Lvi = hat_A{i}' * f{i, q} + aux1 - lambda * f{i, q} + hat_C{i}' * r;
        Cons{i, q} = Lvi <= -e1;
    end
end

Constraints = [Constraints, Cons{:}];

%% (2).2式 # 5-6

Cons = cell(mode_number, divide_number);

for i = 1 : mode_number
    for q = 1 : divide_number
        aux1 = divide_number * (f{i, q + 1} - f{i, q}) / tmin;
        Lvi = hat_A{i}' * f{i, q + 1} + aux1 - lambda * f{i, q + 1} + hat_C{i}' * r;
        Cons{i, q} = Lvi <= -e1; 
    end
end

Constraints = [Constraints, Cons{:}];

%% (2)式 # 7-8

Cons = cell(1, mode_number);

for i = 1 : mode_number 
    Lvi = hat_A{i}' * f{i, divide_number + 1} - lambda * f{i, divide_number + 1} + hat_C{i}' * r;
    Cons{1, i} = Lvi <= -e1; 
end

Constraints = [Constraints, Cons{:}];

%% (3)式 # 9-

Cons = cell(mode_number, divide_number + 1);

for i = 1 : mode_number 
    for q = 1 : divide_number + 1
        Lvi = hat_B{i}' * f{i, q} - gamma * 1 + hat_D{i}' * r;
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
 
% %% 辅助约束
% 
% Cons = cell(mode_number, divide_number + 1);
% 
% for i = 1 : mode_number
%     for j = 1 : divide_number + 1 
%         Lvi = f{i, j};
%         Cons{i, j} = Lvi <= 1000; 
%     end
% end
% 
% Constraints = [Constraints, Cons{:}];


%% 目标函数
 
Objective = gamma;


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
    result.gamma = value(gamma) / xi;
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







