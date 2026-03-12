% =========================================================================
% 双面凸透镜 - 3D 精确光线追迹 
% Author:热带鱼
% Time:2026.3.12
% =========================================================================

clear
clc
close all

%% 1. 定义透镜和环境参数 (你的"全局控制面板")
% 透镜几何参数
R1 = 50;           % 前表面曲率半径 (单位: mm)
R2 = -50;          % 后表面曲率半径 (单位: mm)
d  = 10;           % 透镜中心厚度 (单位: mm)
D  = 30;           % 透镜的有效孔径 (单位: mm)

% 波长与真实玻璃材质定义
wavelength_nm = 587.6; % 入射光波长 (单位: nm)。(486.1nm 蓝光; 587.6nm 黄绿光; 656.3nm 红光)
glass_type = 'N-BK7';  % 透镜玻璃材质 (例如 'N-BK7', 'F2', 'FUSED_SILICA')

% 物理与环境参数
n_env  = 1.0;      % 环境折射率 (空气)
z_start = -20;     % 光线发射的起始 Z 坐标

%光线倾斜入射控制
tilt_angle_deg = 10; % 光线入射倾角 (单位: 度。0 代表正对光轴，大于 0 代表斜射)

% === 1.4 算法控制开关 ===
% use_optimization = true; % true:开启寻优找最小弥散圆 | false:关闭寻优，强制使用理论近轴焦面
use_optimization = false; % true:开启寻优找最小弥散圆 | false:关闭寻优，强制使用理论近轴焦面



%% 1.5根据波长自动计算真实折射率
wv_um = wavelength_nm / 1000;
n_lens = calc_refractive_index(glass_type, wv_um);

fprintf('====================================================\n');
fprintf('当前仿真波长 : %.1f nm\n', wavelength_nm);
fprintf('透镜材质 (%s) 真实折射率 : %.5f\n', glass_type, n_lens);
fprintf('----------------------------------------------------\n');

%% 2. 生成 3D 入射光线束
num_rings = 5;    % 环的数量  
rays_P0 = [];     % 存储所有光线的起点  
rays_V0 = [];     % 存储所有光线的方向  

% 将倾斜角转换为弧度，并计算全局的初始方向向量 (在 Y-Z 平面上向上倾斜)
tilt_rad = tilt_angle_deg * pi / 180;
dir_V = [0, sin(tilt_rad), cos(tilt_rad)];

% 生成中心光线
rays_P0 = [0, 0, z_start];
rays_V0 = dir_V;

% 生成外围环形光线
r_max = D/2 * 0.99; 
for r_idx = 1:num_rings
    radius = r_max * (r_idx / num_rings);%生成一圈又一圈的"环"
    num_pts = 6 * r_idx;  %第 1 环 6 根，第 2 环 12 根，第 3 环 18 根……
    theta = linspace(0, 2*pi, num_pts + 1);
    theta(end) = []; 
    
    %极坐标到直角坐标的转换
    for i = 1:length(theta)
        x = radius * cos(theta(i));
        y = radius * sin(theta(i));
        rays_P0 = [rays_P0; x, y, z_start];%追加到矩阵的末尾。
        rays_V0 = [rays_V0; dir_V]; % 所有光线使用统一的倾斜方向
    end
end
num_total_rays = size(rays_P0, 1);

%% 3. 寻找基准焦平面 (同时计算近轴焦点与边缘焦点) 
% 追迹前表面
C1 = [0, 0, R1];
C2 = [0, 0, d + R2];

% 注意：即使光线斜射，焦平面依然是那个垂直于光轴的感光元件平面
% ------ 3.1 传统的近轴焦点 (使用靠近光轴的光线 y = 0.001) ------
paraxial_P = [0, 0.001, z_start];
paraxial_V = [0, 0, 1]; % 寻找焦平面时，始终用沿光轴的光线寻找

%调用求交函数 intersect_sphere_3d，算出这根探测光线打在前表面的坐标 P1_p
[P1_p, ~] = intersect_sphere_3d(paraxial_P, paraxial_V, C1, R1, 'front');
%算出法向量 N1_p 并确保它迎着光线  
N1_p = (P1_p - C1) / norm(P1_p - C1); if dot(paraxial_V, N1_p) > 0, N1_p = -N1_p; end
% 调用斯涅尔定律函数 snells_law_3d，算出光线进入玻璃后的新方向 V1_p
V1_p = snells_law_3d(paraxial_V, N1_p, n_env, n_lens);

[P2_p, ~] = intersect_sphere_3d(P1_p, V1_p, C2, R2, 'back');
N2_p = (P2_p - C2) / norm(P2_p - C2); if dot(V1_p, N2_p) > 0, N2_p = -N2_p; end
V2_p = snells_law_3d(V1_p, N2_p, n_lens, n_env);

% 计算焦点 Z 坐标 (基于 Y-Z 平面斜率)
z_paraxial_focus = P2_p(3) - P2_p(2) * V2_p(3) / V2_p(2);

% ------ 3.2 边缘焦点 (使用边缘光线高度 y = r_max) ------
marginal_P = [0, r_max, z_start]; % 使用通光孔径的最边缘高度
marginal_V = [0, 0, 1];

[P1_m, hit1_m] = intersect_sphere_3d(marginal_P, marginal_V, C1, R1, 'front');
if hit1_m
    N1_m = (P1_m - C1) / norm(P1_m - C1); if dot(marginal_V, N1_m) > 0, N1_m = -N1_m; end
    V1_m = snells_law_3d(marginal_V, N1_m, n_env, n_lens);
    
    [P2_m, hit2_m] = intersect_sphere_3d(P1_m, V1_m, C2, R2, 'back');
    if hit2_m
        N2_m = (P2_m - C2) / norm(P2_m - C2); if dot(V1_m, N2_m) > 0, N2_m = -N2_m; end
        V2_m = snells_law_3d(V1_m, N2_m, n_lens, n_env);
        
        % 计算边缘光线与光轴的交点 Z 坐标
        z_marginal_focus = P2_m(3) - P2_m(2) * V2_m(3) / V2_m(2);
    end
end

% 打印物理分析数据
fprintf('理论近轴焦点 Z 坐标 (高度 y=0.001) : %.3f mm\n', z_paraxial_focus);
fprintf('实际边缘焦点 Z 坐标 (高度 y=%.3f) : %.3f mm\n', r_max, z_marginal_focus);
fprintf('▶ 系统纵向球差 (LSA)              : %.3f mm\n', z_paraxial_focus - z_marginal_focus);
fprintf('====================================================\n');

% 根据你的需求，寻找平面设置为"近轴焦点"或"边缘光线焦点"，自己调整
z_focus = z_paraxial_focus;  %近轴焦点
% z_focus = z_marginal_focus;  %边缘光线焦点

%% 3.5 最小弥散圆寻优算法 (Best Focus Optimization) 
if use_optimization
% 1. 先追迹一次所有光线，收集它们飞出透镜时的位置(P2)和方向(V2)
valid_P2 = []; valid_V2 = [];
for i = 1:num_total_rays
    P0 = rays_P0(i, :); V0 = rays_V0(i, :);
    [P1, hit1] = intersect_sphere_3d(P0, V0, C1, R1, 'front'); if ~hit1, continue; end
    N1 = (P1 - C1) / norm(P1 - C1); if dot(V0, N1) > 0, N1 = -N1; end
    V1 = snells_law_3d(V0, N1, n_env, n_lens);

    [P2, hit2] = intersect_sphere_3d(P1, V1, C2, R2, 'back'); if ~hit2, continue; end
    N2 = (P2 - C2) / norm(P2 - C2); if dot(V1, N2) > 0, N2 = -N2; end
    V2 = snells_law_3d(V1, N2, n_lens, n_env);

    valid_P2 = [valid_P2; P2]; valid_V2 = [valid_V2; V2];
end

if ~isempty(valid_P2)
    % 2. 定义一个数学目标函数：给定任意 Z 坐标，算出该平面上的 RMS 半径
    % (它调用了写在代码最底部的新增辅助函数 calc_rms_for_z)
    rms_objective_func = @(Z) calc_rms_for_z(Z, valid_P2, valid_V2);
    % 在焦平面附近寻找最佳像面 (现在扫描起点是边缘焦点)
       % 3. 使用 MATLAB 的一维优化器 fminbnd，在近轴焦点附近"扫描"最小值
    % 因为球面像差总是让光线提前汇聚，所以我们往前找15mm，往后找10mm
    [z_best_focus, min_rms] = fminbnd(rms_objective_func, z_focus - 15, z_focus + 10);
    fprintf('▶ 寻优成功！最小弥散圆 Z 坐标   : %.3f mm (最小 RMS: %.4f mm)\n', z_best_focus, min_rms);
    z_focus = z_best_focus;
end
else
    fprintf('▶ 寻优已关闭，强制使用理论近轴焦平面进行作图。\n');
end
%% 4. 执行 3D 光线追迹
%用来记录每一根光线最终打在屏幕上的 X和 Y坐标
spot_x = zeros(num_total_rays, 1);
spot_y = zeros(num_total_rays, 1);

%标记出那些成功走完全程的"有效光线"，画图时把失败的光线剔除掉，当索引用
valid_rays = false(num_total_rays, 1);
figure('Name', sprintf('波长 %.1f nm 3D 透镜追迹', wavelength_nm), 'Position', [100, 100, 1000, 450]);

% --- 子图 1: 3D 光路图 ---
subplot(1, 2, 1);
hold on; grid on; axis equal; view(30, 20);
title(sprintf('3D 光路 (波长 %.1f nm, 倾角 %.1f°)', wavelength_nm, tilt_angle_deg), 'FontSize', 12);
xlabel('X / mm'); ylabel('Y / mm'); zlabel('Z (光轴) / mm');

% 画光轴参考线
plot3([0,0], [0,0], [z_start, z_focus], 'k-.', 'LineWidth', 1);

% 绘制 3D 半透明透镜实体模型
% 1. 生成透镜表面的极坐标网格
[r_grid, theta_grid] = meshgrid(linspace(0, D/2, 30), linspace(0, 2*pi, 40));
x_lens = r_grid .* cos(theta_grid); 
y_lens = r_grid .* sin(theta_grid);

% 2. 计算前后表面的 3D Z坐标 (利用球面方程)
z_lens1 = R1 - sign(R1) * sqrt(R1^2 - r_grid.^2); 
z_lens2 = (d + R2) - sign(R2) * sqrt(R2^2 - r_grid.^2);

% 3. 绘制前后表面 (使用青色 'c'，设置 FaceAlpha 透明度为 0.2)
surf(x_lens, y_lens, z_lens1, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
surf(x_lens, y_lens, z_lens2, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% 4. 绘制透镜外边缘的圆柱侧面 (将前后表面"缝合"成一个实心玻璃块)
theta_edge = linspace(0, 2*pi, 40); 
x_edge = (D/2) * cos(theta_edge); 
y_edge = (D/2) * sin(theta_edge);
z_edge1_val = R1 - sign(R1) * sqrt(R1^2 - (D/2)^2); 
z_edge2_val = (d + R2) - sign(R2) * sqrt(R2^2 - (D/2)^2);

X_edge = [x_edge; x_edge]; 
Y_edge = [y_edge; y_edge]; 
Z_edge = [repmat(z_edge1_val, 1, 40); 
repmat(z_edge2_val, 1, 40)];

surf(X_edge, Y_edge, Z_edge, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

ray_color = 'g-'; if wavelength_nm < 500, ray_color = 'b-'; end; if wavelength_nm > 620, ray_color = 'r-'; end


%光线追迹过程
for i = 1:num_total_rays
    P0 = rays_P0(i, :); V0 = rays_V0(i, :);

    % 第一面
    %调用求交函数寻找它和前表面的交点 P1
    [P1, hit1] = intersect_sphere_3d(P0, V0, C1, R1, 'front'); 
    if ~hit1, continue; end
    % 算出法线 N1 并修正方向
    N1 = (P1 - C1) / norm(P1 - C1); 
    if dot(V0, N1) > 0, N1 = -N1; end
    % 斯涅尔定律计算出进入玻璃内部的新方向 V1
    V1 = snells_law_3d(V0, N1, n_env, n_lens);

    [P2, hit2] = intersect_sphere_3d(P1, V1, C2, R2, 'back'); 
    if ~hit2, continue; end
    N2 = (P2 - C2) / norm(P2 - C2); 
    if dot(V1, N2) > 0, N2 = -N2; end
    V2 = snells_law_3d(V1, N2, n_lens, n_env);
    
    % 投射到焦平面
    t_end = (z_focus - P2(3)) / V2(3); 
    P3 = P2 + t_end * V2;

    spot_x(i) = P3(1); 
    spot_y(i) = P3(2); 
    valid_rays(i) = true;
    
    %绘图
    %起点 到 前表面 (P0 到 P1)
    plot3([P0(1), P1(1)], [P0(2), P1(2)], [P0(3), P1(3)], ray_color, 'LineWidth', 0.5);
    % 前表面 到 后表面 (P1 到 P2)
    plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'c-', 'LineWidth', 0.5);
    % 后表面 到 屏幕 (P2 到 P3)
    plot3([P2(1), P3(1)], [P2(2), P3(2)], [P2(3), P3(3)], ray_color, 'LineWidth', 0.5);
end

% 绘制焦平面边界框
%算出了偏离中心最远的那根光线跑了多远,留出 20% 的视觉边缘空白
plane_size = max(D/2, max(abs(spot_y(valid_rays))) * 1.2);

%红色 (Red),透明度设为 10%,红色的边缘线
fill3([-plane_size, plane_size, plane_size, -plane_size], ...
      [-plane_size, -plane_size, plane_size, plane_size], ...
      [z_focus, z_focus, z_focus, z_focus], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r');

%% 5. 绘制点列图 (Spot Diagram)
subplot(1, 2, 2);
hold on; grid on; axis equal;
title(sprintf('点列图 (波长 %.1f nm)', wavelength_nm), 'FontSize', 12);
xlabel('X 坐标 / mm'); ylabel('Y 坐标 / mm');

valid_x = spot_x(valid_rays); 
valid_y = spot_y(valid_rays);

% 绘制点列图散点
marker_color = 'g'; if wavelength_nm < 500, marker_color = 'b'; end; if wavelength_nm > 620, marker_color = 'r'; end
scatter(valid_x, valid_y, 15, marker_color, 'filled', 'MarkerEdgeColor', 'k');

% 计算光斑的质心 (Centroid)
centroid_x = mean(valid_x); 
centroid_y = mean(valid_y);

% 标出质心位置 (红色的十字)
plot(centroid_x, centroid_y, 'r+', 'MarkerSize', 10, 'LineWidth', 2);

% 设置坐标轴范围 (以质心为中心，方便观察像差形状)
spread_x = max(valid_x) - min(valid_x); 
spread_y = max(valid_y) - min(valid_y);
max_spread = max([spread_x, spread_y, 1e-3]) * 0.8; 
xlim([centroid_x - max_spread, centroid_x + max_spread]); ylim([centroid_y - max_spread, centroid_y + max_spread]);

% 画通过质心的参考线
plot(xlim, [centroid_y, centroid_y], 'k--'); 
plot([centroid_x, centroid_x], ylim, 'k--');

% 升级的 RMS 半径计算 (相对于质心)
% 真实的像斑大小，排除了整体偏移的影响
rms_radius = sqrt(mean((valid_x - centroid_x).^2 + (valid_y - centroid_y).^2));
text(min(xlim) + 0.05*max_spread, max(ylim) - 0.1*max_spread, ...
     sprintf('相对于质心的 RMS 半径: %.4f mm\n整体像偏移: %.2f mm\n测量平面 Z: %.2f mm', rms_radius, centroid_y, z_focus), ...
     'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

%% ========================================================================
% 辅助函数区
% =========================================================================

%折射率计算
% glass_name: 你想查询的玻璃牌号（比如 'N-BK7'）。
% wv_um: 当前入射光的波长，注意这里的单位是微米 (μm)。比如 587.6 nm 的黄绿光，传进来就是 0.5876。
% n: 在这块特定玻璃里的真实折射率。
function n = calc_refractive_index(glass_name, wv_um)
    switch upper(glass_name)
        case 'N-BK7'
            B1 = 1.03961212;  B2 = 0.231792344; B3 = 1.01046945;
            C1 = 0.00600069867; C2 = 0.0200179144; C3 = 103.560653;
        case 'F2'
            B1 = 1.34533359;  B2 = 0.209073176; B3 = 0.937357162;
            C1 = 0.00997743871; C2 = 0.0470450767; C3 = 111.886764;
        case 'FUSED_SILICA'
            B1 = 0.696166300; B2 = 0.407942600; B3 = 0.897479400;
            C1 = 0.00467914826; C2 = 0.0135120631; C3 = 97.9340025;
        otherwise
            n = 1.5; return;
    end
    n2 = 1 + (B1 * wv_um^2)/(wv_um^2 - C1) + (B2 * wv_um^2)/(wv_um^2 - C2) + (B3 * wv_um^2)/(wv_um^2 - C3);
    n = sqrt(n2);
end

%"射线与球面求交算法" (Ray-Sphere Intersection)。
%作用：在三维空间中，已知一条光线的起点和方向，
%以及一个玻璃球面的球心和半径，精确地算出光线打在玻璃上的那个点的 X、Y、Z 坐标。
%%%%%%%%%%%%%%%输入
% P (Point): 光线的起点坐标 $[x, y, z]$。
% V (Vector): 光线飞行的方向向量（必须是长度为 1 的归一化向量）。
% C (Center): 球面的球心坐标 $[x, y, z]$。
% R (Radius): 球面的曲率半径。
% surf_type: 这是一个标签，告诉程序我们要找的是透镜的"前表面" ('front') 还是"后表面" ('back')。
%%%%%%%%%%%%%%%%输出
% P_int: 算出来的精确交点坐标 $[x, y, z]$。
% hit: 一个布尔值 (true/false)，告诉程序这根光线到底有没有射中这个球面。
function [P_int, hit] = intersect_sphere_3d(P, V, C, R, surf_type)
% 向量的点乘（dot）
    OC = P - C;% 从球心指向光线起点的向量
    a = dot(V, V);% 也就是 V 的长度平方，因为 V 归一化了，所以 a 永远等于 1
    b = 2 * dot(V, OC);% 一次项系数
    c = dot(OC, OC) - R^2;% 常数项系数
    
    %如果 $\Delta < 0$，方程无实数解。在物理上意味着光线和球面没有交点
    discriminant = b^2 - 4*a*c;
    if discriminant < 0
        P_int = [0, 0, 0]; hit = false; return;
    end
    
    %计算交点的距离 $t$ (求根公式)
    t1 = (-b - sqrt(discriminant)) / (2*a);
    t2 = (-b + sqrt(discriminant)) / (2*a);
    
    %挑选正确的那个交点 (前表面 vs 后表面)
    if strcmp(surf_type, 'front')
        t = min(t1, t2);% 选距离近的那个点
    else
        t = max(t1, t2);% 选距离远的那个点
    end
    
    %只要飞行距离 $t$ 极小，就认为无效，从而避免浮点数精度问题
    if t < 1e-5
        P_int = [0, 0, 0]; hit = false; return;
    end
    
    %输出最终的三维坐标
    P_int = P + t * V;
    hit = true;
end

%三维向量形式的斯涅尔定律
% V_in: 入射光线的 3D 向量（必须是长度为 1 的单位向量）。
% 
% N: 玻璃表面的 3D 法向量（也必须是单位向量，且必须迎着入射光，即两者夹角大于 90度）。
% 
% n1: 光线当前所在介质的折射率（比如空气，1.0）。
% 
% n2: 光线即将进入介质的折射率（比如玻璃，1.5）。
% 
% V_out: 计算出来的折射（或反射）后的 3D 光线向量。
function V_out = snells_law_3d(V_in, N, n1, n2)
    r = n1 / n2;%准备折射率比值
    cosI = -dot(V_in, N);%计算入射角的余弦值 ($\cos\theta_i$)
    sinT2 = r^2 * (1 - cosI^2);%折射角的正弦平方
    
    %判断是否发生"全反射 (TIR)"
    if sinT2 > 1.0
        V_out = V_in - 2 * dot(V_in, N) * N; % 全反射
    else
        cosT = sqrt(1 - sinT2);
        V_out = r * V_in + (r * cosI - cosT) * N;%3D 向量折射公式
    end
    V_out = V_out / norm(V_out);%归一化
end


%均方根半径 (Root Mean Square Radius)计算
% Z: 优化器当前正在"试探"的感光屏幕的 Z 轴坐标（一个标量数字，比如 50.5）。
% P2: 所有光线从透镜最后一个表面飞出时的起点坐标。它是一个巨大的矩阵（比如 91 * 3），包含了所有 91 根光线的 x, y, z 起点。
% V2: 所有光线飞出时的3D 飞行方向向量。同样是一个 $91 \times 3$ 的大矩阵。
function rms = calc_rms_for_z(Z, P2, V2)
    t = (Z - P2(:,3)) ./ V2(:,3);
    X = P2(:,1) + t .* V2(:,1); Y = P2(:,2) + t .* V2(:,2);
    cx = mean(X); cy = mean(Y);
    rms = sqrt(mean((X - cx).^2 + (Y - cy).^2));
end

%{
Copyright <2026> <热带鱼>
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%}

