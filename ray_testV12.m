% =========================================================================
% 双层凸透镜光线追迹仿真（带MTF分析）
% Author:热带鱼
% Time:2026.3.19
% =========================================================================
clear; clc; close all;

%% 1. 全局参数设置区 (系统核心物理宏观定义) 
wavelength_nm = 587.6;   % 仿真波长 (默认 587.6nm 为氦氖 d 线，光学设计标准参考波长)
z_start = -20;           % 光源发车平面的 Z 轴绝对坐标 (单位: mm)
D = 10;                  % 系统的绝对光阑口径 (Aperture Stop)。在纯透镜系统中，它死死卡在第一面
tilt_angle_deg = 10;      % 光束入射的倾斜角度 (视场角 Field of View)。测试轴上像差设0，偏轴设大角度
num_rings = 5;           % 点列图与光线追迹的同心圆环数 (决定光线的稀密程度)
n_env = 1.0;             % 环境折射率 (默认 1.0 为空气)

% --- 透镜 1 结构参数 ---
Z_center1 = 5;  d1 = 10; R1_1 = 50; R1_2 = -50; glass1 = 'N-BK7'; % 中心Z坐标, 中心厚度, 前表面曲率半径, 后表面曲率半径, 玻璃材质
% --- 透镜 2 结构参数 ---
Z_center2 = 25; d2 = 10; R2_1 = 60; R2_2 = -60; glass2 = 'F2';    

use_optimization = false; % 永久关闭寻优，死守理论近轴焦平面，以此来直观暴露大视场下的场曲和离焦像差
N_grid = 256;             % 波动光学 FFT 采样的网格分辨率 (256x256，决定波前图和 PSF 的清晰度)

%% 1.5 物理引擎：根据波长自动计算真实折射率 (色散模型) 
% 将纳米转换为微米，代入底部的 Sellmeier 色散公式，获取精准的材料折射率
wv_um = wavelength_nm / 1000;
n_lens1 = calc_refractive_index(glass1, wv_um);
n_lens2 = calc_refractive_index(glass2, wv_um);

%% 2. 内部计算初始化 (构建 3D 几何模型坐标)
% 计算透镜前后表面的顶点 (Vertex) 坐标
z_v1_front = Z_center1 - d1/2; z_v1_back  = Z_center1 + d1/2;
C1 = [0, 0, z_v1_front + R1_1]; C2 = [0, 0, z_v1_back  + R1_2]; % C1, C2 为透镜1前后表面的球心坐标

z_v2_front = Z_center2 - d2/2; z_v2_back  = Z_center2 + d2/2;
C3 = [0, 0, z_v2_front + R2_1]; C4 = [0, 0, z_v2_back  + R2_2]; % C3, C4 为透镜2前后表面的球心坐标

% 【核心对准逻辑】：计算大视场斜射时的发车偏移量
% 为了确保倾斜光束在飞行 z_start 到 z_v1_front 的距离后，依然能**正中**第一片透镜的中心(0,0)，
% 必须根据倾角和飞行距离，反向推算出发车点需要的 X/Y 偏移补偿。
tilt_rad = tilt_angle_deg * pi / 180;
dir_V = [0, sin(tilt_rad), cos(tilt_rad)]; % 初始光束的统一方向向量
delta_z = z_v1_front - z_start;            % 飞行距离
shift_x = delta_z * dir_V(1) / dir_V(3);   % X 方向漂移量
shift_y = delta_z * dir_V(2) / dir_V(3);   % Y 方向漂移量

% 恒定入瞳策略：光束永远是 D/2 大小，不再收缩！(0.95 是为了防止刚好压在边界导致插值计算崩溃)
r_max = D/2 * 0.95; 

%% 3. 生成 3D 入射光线束 (点列图专用几何光束)
rays_P0 = []; rays_V0 = [];       
rays_P0 = [-shift_x, -shift_y, z_start]; rays_V0 = dir_V; % 注入中心主光线

% 循环生成同心圆环阵列的光线
for r_idx = 1:num_rings
    radius = r_max * (r_idx / num_rings);
    num_pts = 6 * r_idx; % 外圈光线更密，保持面积能量均匀
    theta = linspace(0, 2*pi, num_pts + 1); theta(end) = []; 
    for i = 1:length(theta)
        x_pupil = radius * cos(theta(i)); y_pupil = radius * sin(theta(i));
        % 发射坐标均减去 shift_x/y，确保光束抵达第一面时是一个完美的圆
        rays_P0 = [rays_P0; x_pupil - shift_x, y_pupil - shift_y, z_start]; 
        rays_V0 = [rays_V0; dir_V]; 
    end
end
num_total_rays = size(rays_P0, 1);

%% 4. 寻找精确的"近轴焦平面" 与 "纵向球差(LSA)"
% 近轴焦点定位：发射一根高度仅为 1 微米 (0.001mm) 的极细光线，它不会产生球差，代表最理想的理论高斯焦点
paraxial_P = [0, 0.001, z_start]; paraxial_V = [0, 0, 1]; 
[P1_p, ~] = intersect_sphere_3d(paraxial_P, paraxial_V, C1, R1_1, 'front');
N1_p = (P1_p - C1) / norm(P1_p - C1); if dot(paraxial_V, N1_p) > 0, N1_p = -N1_p; end; V1_p = snells_law_3d(paraxial_V, N1_p, n_env, n_lens1);
[P2_p, ~] = intersect_sphere_3d(P1_p, V1_p, C2, R1_2, 'back');
N2_p = (P2_p - C2) / norm(P2_p - C2); if dot(V1_p, N2_p) > 0, N2_p = -N2_p; end; V2_p = snells_law_3d(V1_p, N2_p, n_lens1, n_env);
[P3_p, ~] = intersect_sphere_3d(P2_p, V2_p, C3, R2_1, 'front');
N3_p = (P3_p - C3) / norm(P3_p - C3); if dot(V2_p, N3_p) > 0, N3_p = -N3_p; end; V3_p = snells_law_3d(V2_p, N3_p, n_env, n_lens2);
[P4_p, ~] = intersect_sphere_3d(P3_p, V3_p, C4, R2_2, 'back');
N4_p = (P4_p - C4) / norm(P4_p - C4); if dot(V3_p, N4_p) > 0, N4_p = -N4_p; end; V4_p = snells_law_3d(V3_p, N4_p, n_lens2, n_env);
z_focus = P4_p(3) - P4_p(2) * V4_p(3) / V4_p(2); % 相似三角形外推求交 Z 轴

% 附加分析：追踪一根轴上最边缘光线 (高度 D/2)，计算纵向球差 (LSA)
marg_P = [0, D/2 * 0.99, z_start]; marg_V = [0, 0, 1];
[P1_m, hit_m] = intersect_sphere_3d(marg_P, marg_V, C1, R1_1, 'front');
if hit_m
    N1_m = (P1_m - C1) / norm(P1_m - C1); if dot(marg_V, N1_m) > 0, N1_m = -N1_m; end; V1_m = snells_law_3d(marg_V, N1_m, n_env, n_lens1);
    [P2_m, ~] = intersect_sphere_3d(P1_m, V1_m, C2, R1_2, 'back');
    N2_m = (P2_m - C2) / norm(P2_m - C2); if dot(V1_m, N2_m) > 0, N2_m = -N2_m; end; V2_m = snells_law_3d(V1_m, N2_m, n_lens1, n_env);
    [P3_m, ~] = intersect_sphere_3d(P2_m, V2_m, C3, R2_1, 'front');
    N3_m = (P3_m - C3) / norm(P3_m - C3); if dot(V2_m, N3_m) > 0, N3_m = -N3_m; end; V3_m = snells_law_3d(V2_m, N3_m, n_env, n_lens2);
    [P4_m, ~] = intersect_sphere_3d(P3_m, V3_m, C4, R2_2, 'back');
    N4_m = (P4_m - C4) / norm(P4_m - C4); if dot(V3_m, N4_m) > 0, N4_m = -N4_m; end; V4_m = snells_law_3d(V3_m, N4_m, n_lens2, n_env);
    z_marginal = P4_m(3) - P4_m(2) * V4_m(3) / V4_m(2);
    LSA = z_marginal - z_focus; % 边缘光线焦点与近轴焦点的差值即为球差
else
    LSA = NaN;
end

%% =========================================================================
% 6. 绘图 Figure 1：几何光学分析 (3D 光路 + Zemax同款自适应透镜渲染)
% =========================================================================
spot_x = zeros(num_total_rays, 1); spot_y = zeros(num_total_rays, 1); valid_rays = false(num_total_rays, 1);
figure('Name', sprintf('波长 %.1f nm 几何光学分析', wavelength_nm), 'Position', [50, 100, 1000, 450]);
subplot(1, 2, 1); hold on; grid on; axis equal; view(30, 20);
title(sprintf('双透镜 3D 光路 (Tilt: %.1f°)', tilt_angle_deg), 'FontSize', 12);
xlabel('X / mm'); ylabel('Y / mm'); zlabel('Z / mm');
plot3([0,0], [0,0], [z_start, z_focus], 'k-.', 'LineWidth', 1); 


% 在正式画图前，先闭着眼睛把所有光线跑一遍，记录它们在每片透镜上留下的最大扩散足迹。
D_relax = D * 4; % 物理计算时放宽限制，确保探测到最极限的光束尺寸
r_max_L1 = D/2;  % 第一面是光阑面，其机械边界至少等于光阑口径
r_max_L2 = 0.1;  % 初始化透镜2的探测下限
for i = 1:num_total_rays
    P0 = rays_P0(i, :); V0 = rays_V0(i, :);
    [P1, hit1] = intersect_sphere_3d(P0, V0, C1, R1_1, 'front'); 
    if ~hit1 || (P1(1)^2 + P1(2)^2 > (D/2)^2), continue; end % 面1是系统光阑，在此硬截断
    N1 = (P1 - C1) / norm(P1 - C1); if dot(V0, N1) > 0, N1 = -N1; end; V1 = snells_law_3d(V0, N1, n_env, n_lens1);
    
    [P2, hit2] = intersect_sphere_3d(P1, V1, C2, R1_2, 'back'); 
    if ~hit2 || (P2(1)^2 + P2(2)^2 > (D_relax/2)^2), continue; end
    r_max_L1 = max(r_max_L1, norm(P2(1:2))); % 探测透镜1后表面的最大极限径向坐标
    N2 = (P2 - C2) / norm(P2 - C2); if dot(V1, N2) > 0, N2 = -N2; end; V2 = snells_law_3d(V1, N2, n_lens1, n_env);
    
    [P3, hit3] = intersect_sphere_3d(P2, V2, C3, R2_1, 'front'); 
    if ~hit3 || (P3(1)^2 + P3(2)^2 > (D_relax/2)^2), continue; end
    r_max_L2 = max(r_max_L2, norm(P3(1:2))); % 探测透镜2前表面的最大极限
    N3 = (P3 - C3) / norm(P3 - C3); if dot(V2, N3) > 0, N3 = -N3; end; V3 = snells_law_3d(V2, N3, n_env, n_lens2);
    
    [P4, hit4] = intersect_sphere_3d(P3, V3, C4, R2_2, 'back'); 
    if ~hit4 || (P4(1)^2 + P4(2)^2 > (D_relax/2)^2), continue; end
    r_max_L2 = max(r_max_L2, norm(P4(1:2))); % 探测透镜2后表面的最大极限
end

% 为透镜实体增加 2% 的机械框余量，避免光线紧贴边界，提升视觉美感
r_draw_L1 = r_max_L1 * 1.02; 
r_draw_L2 = r_max_L2 * 1.02;

% 实体渲染 1：精准包裹透镜 1 (包含前后球面与侧向圆柱面)
[rg1, tg1] = meshgrid(linspace(0, r_draw_L1, 30), linspace(0, 2*pi, 40));
xl1 = rg1 .* cos(tg1); yl1 = rg1 .* sin(tg1);
z1_f = z_v1_front + R1_1 - sign(R1_1)*sqrt(R1_1^2 - rg1.^2); % 矢高公式计算球面 Z 坐标
z1_b = z_v1_back  + R1_2 - sign(R1_2)*sqrt(R1_2^2 - rg1.^2);
surf(xl1, yl1, z1_f, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
surf(xl1, yl1, z1_b, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% 透镜 1 紧凑侧边 (连接前后圆周)
xe1 = r_draw_L1 * cos(tg1(:,1)); ye1 = r_draw_L1 * sin(tg1(:,1));
z1_fe = z_v1_front + R1_1 - sign(R1_1)*sqrt(R1_1^2 - r_draw_L1^2);
z1_be = z_v1_back  + R1_2 - sign(R1_2)*sqrt(R1_2^2 - r_draw_L1^2);
surf([xe1, xe1], [ye1, ye1], [repmat(z1_fe, 40, 1), repmat(z1_be, 40, 1)], 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.15);

% 渲染独立的金属系统光阑 (Aperture Stop)
% 只有当光线漂移导致透镜实体半径大于 D/2 时，才画出这个遮挡多余光线的黑环
if r_draw_L1 > D/2
    [rg_stop, tg_stop] = meshgrid(linspace(D/2, r_draw_L1, 10), linspace(0, 2*pi, 40));
    xs = rg_stop .* cos(tg_stop); ys = rg_stop .* sin(tg_stop);
    zs = z_v1_front + R1_1 - sign(R1_1)*sqrt(R1_1^2 - rg_stop.^2) + 0.05; % 紧贴前表面悬浮
    surf(xs, ys, zs, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.6); 
end

% 实体渲染 2：精准包裹透镜 2 (同理)
[rg2, tg2] = meshgrid(linspace(0, r_draw_L2, 30), linspace(0, 2*pi, 40));
xl2 = rg2 .* cos(tg2); yl2 = rg2 .* sin(tg2);
z2_f = z_v2_front + R2_1 - sign(R2_1)*sqrt(R2_1^2 - rg2.^2);
z2_b = z_v2_back  + R2_2 - sign(R2_2)*sqrt(R2_2^2 - rg2.^2);
surf(xl2, yl2, z2_f, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
surf(xl2, yl2, z2_b, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.15);
xe2 = r_draw_L2 * cos(tg2(:,1)); ye2 = r_draw_L2 * sin(tg2(:,1));
z2_fe = z_v2_front + R2_1 - sign(R2_1)*sqrt(R2_1^2 - r_draw_L2^2);
z2_be = z_v2_back  + R2_2 - sign(R2_2)*sqrt(R2_2^2 - r_draw_L2^2);
surf([xe2, xe2], [ye2, ye2], [repmat(z2_fe, 40, 1), repmat(z2_be, 40, 1)], 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.15);

% ================= 正式画光路 =================
for i = 1:num_total_rays
    P0 = rays_P0(i, :); V0 = rays_V0(i, :);
    
    % 【光阑面】：严格卡死进光量 D/2
    [P1, hit1] = intersect_sphere_3d(P0, V0, C1, R1_1, 'front'); 
    if ~hit1 || (P1(1)^2 + P1(2)^2 > (D/2)^2), continue; end
    N1 = (P1 - C1) / norm(P1 - C1); if dot(V0, N1) > 0, N1 = -N1; end; V1 = snells_law_3d(V0, N1, n_env, n_lens1);
    
    % 【放宽面】：面2~面4 纯负责折射，使用 D_relax 放宽物理约束，避免渐晕误杀
    [P2, hit2] = intersect_sphere_3d(P1, V1, C2, R1_2, 'back'); 
    if ~hit2 || (P2(1)^2 + P2(2)^2 > (D_relax/2)^2), continue; end
    N2 = (P2 - C2) / norm(P2 - C2); if dot(V1, N2) > 0, N2 = -N2; end; V2 = snells_law_3d(V1, N2, n_lens1, n_env);
    
    [P3, hit3] = intersect_sphere_3d(P2, V2, C3, R2_1, 'front'); 
    if ~hit3 || (P3(1)^2 + P3(2)^2 > (D_relax/2)^2), continue; end
    N3 = (P3 - C3) / norm(P3 - C3); if dot(V2, N3) > 0, N3 = -N3; end; V3 = snells_law_3d(V2, N3, n_env, n_lens2);
    
    [P4, hit4] = intersect_sphere_3d(P3, V3, C4, R2_2, 'back'); 
    if ~hit4 || (P4(1)^2 + P4(2)^2 > (D_relax/2)^2), continue; end
    N4 = (P4 - C4) / norm(P4 - C4); if dot(V3, N4) > 0, N4 = -N4; end; V4 = snells_law_3d(V3, N4, n_lens2, n_env);
    
    t_end = (z_focus - P4(3)) / V4(3); P5 = P4 + t_end * V4; % 飞向接收靶面
    spot_x(i) = P5(1); spot_y(i) = P5(2); valid_rays(i) = true; % 记录有效像点
    
    % 绘制各段光线，展现穿过界面的折射路径
    plot3([P0(1), P1(1)], [P0(2), P1(2)], [P0(3), P1(3)], 'g-', 'LineWidth', 0.5); 
    plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'c-', 'LineWidth', 0.5); 
    plot3([P2(1), P3(1)], [P2(2), P3(2)], [P2(3), P3(3)], 'g-', 'LineWidth', 0.5); 
    plot3([P3(1), P4(1)], [P3(2), P4(2)], [P3(3), P4(3)], 'm-', 'LineWidth', 0.5); 
    plot3([P4(1), P5(1)], [P4(2), P5(2)], [P4(3), P5(3)], 'g-', 'LineWidth', 0.5); 
end
% 画出接收面 (红色靶板)
plane_size = max(D/2, max(abs(spot_y(valid_rays))) * 1.2);
fill3([-plane_size, plane_size, plane_size, -plane_size], [-plane_size, -plane_size, plane_size, plane_size], [z_focus, z_focus, z_focus, z_focus], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'r');

% ================= 绘制点列图子图 =================
subplot(1, 2, 2); hold on; grid on; axis equal;
title('点列图 (Spot Diagram)', 'FontSize', 12);
xlabel('X 坐标 / mm'); ylabel('Y 坐标 / mm');
valid_x = spot_x(valid_rays); valid_y = spot_y(valid_rays);
if isempty(valid_x), text(0, 0, '所有光线均已阵亡', 'Color', 'r', 'HorizontalAlignment', 'center'); return; end
scatter(valid_x, valid_y, 15, 'g', 'filled', 'MarkerEdgeColor', 'k');
centroid_x = mean(valid_x); centroid_y = mean(valid_y); % 寻找能量质心
plot(centroid_x, centroid_y, 'r+', 'MarkerSize', 10, 'LineWidth', 2);

% 【恢复数据仪表盘】：计算高级光学评价指标
rms_radius = sqrt(mean((valid_x - centroid_x).^2 + (valid_y - centroid_y).^2)); % 均方根弥散圆半径
spread_x = max(valid_x) - min(valid_x); spread_y = max(valid_y) - min(valid_y);
max_spread = max([spread_x, spread_y, 1e-3]) * 0.8; 
xlim([centroid_x - max_spread, centroid_x + max_spread]); ylim([centroid_y - max_spread, centroid_y + max_spread]);

% 在图表左上角挂载高级数据面板 (包含 Z坐标, RMS, 像移, LSA球差)
x_lims = xlim; y_lims = ylim;
info_text = sprintf('测量焦平面 Z : %.4f mm\nRMS 半径 : %.4f mm\n像偏移 Y : %.4f mm\n纵向球差 LSA : %.4f mm', ...
                    z_focus, rms_radius, centroid_y, LSA);
text(x_lims(1) + 0.05*max_spread, y_lims(2) - 0.1*max_spread, info_text, ...
     'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'k');

% 并在控制台同步打印输出报告
fprintf('▶ 测量平面 Z 坐标 (近轴焦点) : %.5f mm\n', z_focus);
fprintf('▶ 纵向球差 (LSA)             : %.5f mm\n', LSA);
fprintf('▶ RMS 光斑半径               : %.5f mm\n', rms_radius);
fprintf('▶ 质心像偏移 Y               : %.5f mm\n', centroid_y);
fprintf('====================================================\n');

%% =========================================================================
% 8. 绘图 Figure 2：物理波动光学评价 (Wavefront, PSF & FFT MTF) 
% =========================================================================
figure('Name', sprintf('波长 %.1f nm 物理波动光学分析', wavelength_nm), 'Position', [50, 150, 1500, 450]);

% 网格必须死死锚定入瞳口径 D，提取纯净的原始相位
x_vec = linspace(-D/2, D/2, N_grid); y_vec = linspace(-D/2, D/2, N_grid);
[X_pupil, Y_pupil] = meshgrid(x_vec, y_vec);
dx_pupil = x_vec(2) - x_vec(1);
Pupil_Function = zeros(N_grid, N_grid); 
Wavefront_Map = NaN(N_grid, N_grid); 

% === 计算中心主光线的绝对参考光程 (标尺) ===
P0_c = [-shift_x, -shift_y, z_start]; V0_c = dir_V;
[P1_c,~,t1_c] = intersect_sphere_3d(P0_c, V0_c, C1, R1_1, 'front'); N1_c = (P1_c - C1)/norm(P1_c - C1); if dot(V0_c, N1_c)>0, N1_c = -N1_c; end; V1_c = snells_law_3d(V0_c, N1_c, n_env, n_lens1);
[P2_c,~,t2_c] = intersect_sphere_3d(P1_c, V1_c, C2, R1_2, 'back');  N2_c = (P2_c - C2)/norm(P2_c - C2); if dot(V1_c, N2_c)>0, N2_c = -N2_c; end; V2_c = snells_law_3d(V1_c, N2_c, n_lens1, n_env);
[P3_c,~,t3_c] = intersect_sphere_3d(P2_c, V2_c, C3, R2_1, 'front'); N3_c = (P3_c - C3)/norm(P3_c - C3); if dot(V2_c, N3_c)>0, N3_c = -N3_c; end; V3_c = snells_law_3d(V2_c, N3_c, n_env, n_lens2);
[P4_c,~,t4_c] = intersect_sphere_3d(P3_c, V3_c, C4, R2_2, 'back');  N4_c = (P4_c - C4)/norm(P4_c - C4); if dot(V3_c, N4_c)>0, N4_c = -N4_c; end; V4_c = snells_law_3d(V3_c, N4_c, n_lens2, n_env);
t5_c = (z_focus - P4_c(3)) / V4_c(3); P5_c = P4_c + t5_c * V4_c;
F_ref = P5_c; % 建立理想参考球面波的球心
% 计算主光线的总光程，并将其作为减去相位常量的基准
OPL_chief_screen = n_env*t1_c + n_lens1*t2_c + n_env*t3_c + n_lens2*t4_c + n_env*t5_c; 
OPL_chief_ref = OPL_chief_screen + n_env * dot(V4_c, F_ref - P5_c);
lambda_mm = wavelength_nm * 1e-6;

% 【核心物理雷达】：追踪真实的有效出瞳数值孔径 (NA)
NA_max = 0;      % 存储能够成功穿越系统边缘的真实像方角
r_ent_max = 0;   % 记录该极限光线在入瞳上的起始坐标

% ================== 暴力网格追迹累积光程差 ==================
for ii = 1:N_grid
    for jj = 1:N_grid
        x0_p = X_pupil(ii, jj); y0_p = Y_pupil(ii, jj);
        if x0_p^2 + y0_p^2 > (D/2)^2, continue; end % 切除光阑外的光线
        
        P0 = [x0_p - shift_x, y0_p - shift_y, z_start]; V0 = dir_V;
        
        [P1, hit1, t1] = intersect_sphere_3d(P0, V0, C1, R1_1, 'front'); 
        if ~hit1 || (P1(1)^2 + P1(2)^2 > (D/2)^2), continue; end
        N1 = (P1 - C1)/norm(P1 - C1); if dot(V0, N1)>0, N1 = -N1; end; V1 = snells_law_3d(V0, N1, n_env, n_lens1);
        
        [P2, hit2, t2] = intersect_sphere_3d(P1, V1, C2, R1_2, 'back'); 
        if ~hit2 || (P2(1)^2 + P2(2)^2 > (D_relax/2)^2), continue; end
        N2 = (P2 - C2)/norm(P2 - C2); if dot(V1, N2)>0, N2 = -N2; end; V2 = snells_law_3d(V1, N2, n_lens1, n_env);
        
        [P3, hit3, t3] = intersect_sphere_3d(P2, V2, C3, R2_1, 'front'); 
        if ~hit3 || (P3(1)^2 + P3(2)^2 > (D_relax/2)^2), continue; end
        N3 = (P3 - C3)/norm(P3 - C3); if dot(V2, N3)>0, N3 = -N3; end; V3 = snells_law_3d(V2, N3, n_env, n_lens2);
        
        [P4, hit4, t4] = intersect_sphere_3d(P3, V3, C4, R2_2, 'back'); 
        if ~hit4 || (P4(1)^2 + P4(2)^2 > (D_relax/2)^2), continue; end
        N4 = (P4 - C4)/norm(P4 - C4); if dot(V3, N4)>0, N4 = -N4; end; V4 = snells_law_3d(V3, N4, n_lens2, n_env);
        
        t5 = (z_focus - P4(3)) / V4(3); P5 = P4 + t5 * V4;
        
        % 【雷达工作】：记录到达焦平面的该光线与中心光线的真实夹角 NA
        r_ent = sqrt(x0_p^2 + y0_p^2);
        cos_theta = dot(V4, V4_c);
        sin_theta = sqrt(max(0, 1 - cos_theta^2));
        if sin_theta > NA_max
            NA_max = sin_theta; % 更新最大数值孔径
            r_ent_max = r_ent;  % 记录发生该 NA 的对应入瞳物理坐标
        end
        
        % 倾斜面初始光程补偿，确保提取平滑的波前
        initial_OPL = n_env * (x0_p * V0(1) + y0_p * V0(2)); 
        OPL_ray_screen = initial_OPL + n_env*t1 + n_lens1*t2 + n_env*t3 + n_lens2*t4 + n_env*t5;
        OPL_ray_ref = OPL_ray_screen + n_env * dot(V4, F_ref - P5); % 投影回参考球面
        
        OPD = OPL_ray_ref - OPL_chief_ref; % 求解相对光程差
        W_waves = OPD / lambda_mm;         % 转换为波长数量 (Waves)
        
        Wavefront_Map(ii, jj) = W_waves; 
        Pupil_Function(ii, jj) = exp(1i * 2 * pi * W_waves); % 构建复振幅光瞳矩阵
    end
end

% ================== 核心 FFT 计算 (求解 PSF 与 MTF) ==================
N_pad = N_grid * 2; % 补零以提升频率域采样密度 (Zero-padding)
P_pad = zeros(N_pad, N_pad);
start_idx = N_pad/2 - N_grid/2 + 1; end_idx = N_pad/2 + N_grid/2;
P_pad(start_idx:end_idx, start_idx:end_idx) = Pupil_Function;
% 傅里叶变换求光强分布 (PSF)
PSF = abs(fftshift(fft2(ifftshift(P_pad)))).^2;
PSF = PSF / sum(PSF(:)); 
% 对 PSF 进行逆傅里叶变换求光学传递函数 (OTF) 与模值 (MTF)
OTF = fftshift(ifft2(ifftshift(PSF)));
MTF = abs(OTF);
MTF = MTF / max(MTF(:));
R_ref = norm(F_ref - P4_c); 

% 【极其关键的坐标域放缩技术】：用真实探测到的出瞳 NA_max 缩放频率系
if r_ent_max == 0, r_ent_max = D/2; NA_max = (D/2)/R_ref; end % 防止空跑导致的除零报错
% 将入瞳物理空间映射为真实出瞳的角域 (Direction Cosine)
Delta_alpha = dx_pupil * (NA_max / r_ent_max);
% 生成符合物理法则的空间频率横坐标 (X 轴标尺)
dx_psf = lambda_mm / (N_pad * Delta_alpha); 
psf_axis_um = (-(N_pad/2):(N_pad/2 - 1)) * dx_psf * 1000; 
df = Delta_alpha / lambda_mm; % FFT 单个像素代表的 lp/mm
center_idx = N_pad/2 + 1;
f_axis = (0:(N_pad/2 - 1)) * df; % 最终的 MTF 物理刻度尺

% === 子图 1：渲染 3D 坑洼波前 ===
subplot(1, 3, 1); hold on; grid on;
title(sprintf('波前相位图 (Tilt: %.1f°)', tilt_angle_deg), 'FontSize', 12);
xlabel('入瞳 X / mm'); ylabel('入瞳 Y / mm'); zlabel('OPD (Waves)');
surf(X_pupil, Y_pupil, Wavefront_Map, 'EdgeColor', 'none');
shading interp; colormap(gca, 'jet'); colorbar; axis square; view(45, 45); 

% === 子图 2：渲染带有衍射环的点扩散函数 (PSF) ===
subplot(1, 3, 2); hold on; grid on;
title('高清点扩散函数 (PSF)', 'FontSize', 12);
xlabel('X (\mum)'); ylabel('Y (\mum)');
PSF_vis = (PSF / max(PSF(:))).^0.6; % 施加 Gamma 0.6 提亮暗弱的衍射环和彗星尾巴
zoom_N = 128; idx_zoom = (center_idx - zoom_N/2) : (center_idx + zoom_N/2 - 1);
imagesc(psf_axis_um(idx_zoom), psf_axis_um(idx_zoom), PSF_vis(idx_zoom, idx_zoom));
colormap(gca, 'hot'); colorbar; axis tight; axis square; set(gca, 'YDir', 'normal'); 

% === 子图 3：渲染 MTF 与通用衍射极限曲线 ===
subplot(1, 3, 3); hold on; grid on;
title('波动光学 FFT MTF', 'FontSize', 12);
xlabel('空间频率 (lp/mm)'); ylabel('MTF 模值');
MTF_S = MTF(center_idx, center_idx:end); MTF_T = MTF(center_idx:end, center_idx); 
plot(f_axis, MTF_S, 'b-', 'LineWidth', 2, 'DisplayName', 'Sagittal (X)');
plot(f_axis, MTF_T, 'r-', 'LineWidth', 2, 'DisplayName', 'Tangential (Y)');

% 【自相关原理计算物理】：利用真实 NA 计算阿贝理论截止频率
fc = 2 * NA_max / lambda_mm; 
% 强行抹除所有相位误差 (OPD=0)，仅保留透光孔径形状，算出该口径最完美的状态
Pupil_Ideal = abs(Pupil_Function); 
P_pad_ideal = zeros(N_pad, N_pad); P_pad_ideal(start_idx:end_idx, start_idx:end_idx) = Pupil_Ideal;
PSF_ideal = abs(fftshift(fft2(ifftshift(P_pad_ideal)))).^2; PSF_ideal = PSF_ideal / sum(PSF_ideal(:)); 
OTF_ideal = fftshift(ifft2(ifftshift(PSF_ideal))); MTF_ideal = abs(OTF_ideal); MTF_ideal = MTF_ideal / max(MTF_ideal(:));
MTF_diff_universal = MTF_ideal(center_idx, center_idx:end); 
plot(f_axis, MTF_diff_universal, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Diffraction Limit');
legend('Location', 'northeast'); 
xlim([0, min(1000, fc)]); 
ylim([0, 1.05]);

%% =========================================================================
% 底层辅助数学函数区 (折射率、球面交点、斯涅尔定律、RMS算子)
% =========================================================================
function n = calc_refractive_index(glass_name, wv_um)
    % 经典的 Sellmeier 色散公式，根据材料和波长计算绝对折射率
    switch upper(glass_name)
        case 'N-BK7', B1 = 1.03961212; B2 = 0.231792344; B3 = 1.01046945; C1 = 0.00600069867; C2 = 0.0200179144; C3 = 103.560653;
        case 'F2', B1 = 1.34533359; B2 = 0.209073176; B3 = 0.937357162; C1 = 0.00997743871; C2 = 0.0470450767; C3 = 111.886764;
        otherwise, n = 1.5; return;
    end
    n = sqrt(1 + (B1*wv_um^2)/(wv_um^2-C1) + (B2*wv_um^2)/(wv_um^2-C2) + (B3*wv_um^2)/(wv_um^2-C3));
end

function [P_int, hit, t] = intersect_sphere_3d(P, V, C, R, surf_type)
    % 基于二次方程判别式求光线与三维球面的精确数学交点
    OC = P - C; a = dot(V, V); b = 2 * dot(V, OC); c = dot(OC, OC) - R^2;
    discriminant = b^2 - 4*a*c; 
    if discriminant < 0, P_int = [0, 0, 0]; hit = false; t = 0; return; end
    t1 = (-b - sqrt(discriminant))/(2*a); t2 = (-b + sqrt(discriminant))/(2*a);
    if strcmp(surf_type, 'front'), t = min(t1, t2); else, t = max(t1, t2); end
    if t < 1e-6, P_int = [0, 0, 0]; hit = false; t = 0; return; end
    P_int = P + t * V; hit = true;
end

function V_out = snells_law_3d(V_in, N, n1, n2)
    % 三维向量版斯涅尔折射定律计算出射光线向量
    r = n1 / n2; cosI = -dot(V_in, N); sinT2 = r^2 * (1 - cosI^2);
    if sinT2 > 1.0, V_out = V_in - 2 * dot(V_in, N) * N; % 触发全反射现象
    else, cosT = sqrt(1 - sinT2); V_out = r * V_in + (r * cosI - cosT) * N; end
    V_out = V_out / norm(V_out); 
end

function rms = calc_rms_for_z(Z, P2, V2)
    % 给定评估屏幕 Z，计算投影点列图的能量集中度 (均方根半径 RMS)
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




