clc
clear all

% 物理常量
Constant_rho = 1000;
Constant_g = 9.8;

% Wigley船型主尺度信息
Wigley_L = 2;
Wigley_B = 1;
Wigley_D = 0.12504;
Function_y = @(x,z) 0.5 * Wigley_B * (1 - 4 * x.^2 / Wigley_L^2).* (1 - z.^2 / Wigley_D^2);
Wigley_V_Real_half = integral2(Function_y, -0.5 * Wigley_L,0.5 * Wigley_L,-Wigley_D,0);
Wigley_V_Real = 2 * Wigley_V_Real_half;
Function_l = @(x,z) sqrt(1 + ((4 * Wigley_B / Wigley_L.^2).* (z.^2 / Wigley_D.^2 - 1).* x).^2);
Wigley_S_Real_half = integral2(Function_l, -0.5 * Wigley_L,0.5 * Wigley_L,-Wigley_D,0);
Wigley_S_Real = 2 * Wigley_S_Real_half;

% 离散程度设置，Iteration_Wigley_L为沿船长方向等分数量，Iteration_Wigley_D为沿型深方向等分数量
Iteration_Wigley_L = 40;
Iteration_Wigley_D = 10;
delta_Wigley_L = Wigley_L / Iteration_Wigley_L;
delta_Wigley_D = Wigley_D / Iteration_Wigley_D;

% 分别构造XYZ坐标值的矩阵
X_all = linspace(Wigley_L / 2, -Wigley_L / 2,Iteration_Wigley_L + 1);
Z_all = linspace(0, - Wigley_D,Iteration_Wigley_D + 1);
Y_all = zeros(Iteration_Wigley_D + 1,Iteration_Wigley_L + 1);
for Iteration_i = (1:Iteration_Wigley_L + 1)
    for Iteration_j = (1:Iteration_Wigley_D + 1)
        Y_all(Iteration_j,Iteration_i) = 0.5 * Wigley_B * (1 - (X_all(1,Iteration_i) / (Wigley_L * 0.5))^2) * (1 - (Z_all(1,Iteration_j) / (Wigley_D))^2);
    end
end

% 求离散后湿表面积和排水体积
S_Discrete_half = 0;
V_Discrete_half = 0;
for Iteration_m = (1:Iteration_Wigley_L)
    for Iteration_n = (1:Iteration_Wigley_D)
        % 求各点坐标值及边长
        Point_A_x = X_all(1,Iteration_m);
        Point_A_y = Y_all(Iteration_n,Iteration_m);
        Point_A_z = Z_all(1,Iteration_n);
        Point_B_x = X_all(1,Iteration_m);
        Point_B_y = Y_all(Iteration_n + 1,Iteration_m);
        Point_B_z = Z_all(1,Iteration_n + 1);
        Point_C_x = X_all(1,Iteration_m + 1);
        Point_C_y = Y_all(Iteration_n,Iteration_m + 1);
        Point_C_z = Z_all(1,Iteration_n + 1);
        Point_D_x = X_all(1,Iteration_m + 1);
        Point_D_y = Y_all(Iteration_n + 1,Iteration_m + 1);
        Point_D_z = Z_all(1,Iteration_n + 1);
        deltaL_AB = sqrt((Point_A_x - Point_B_x)^2 + (Point_A_y - Point_B_y)^2 + (Point_A_z - Point_B_z)^2);
        deltaL_AC = sqrt((Point_A_x - Point_C_x)^2 + (Point_A_y - Point_C_y)^2 + (Point_A_z - Point_C_z)^2);
        deltaL_BD = sqrt((Point_B_x - Point_D_x)^2 + (Point_B_y - Point_D_y)^2 + (Point_B_z - Point_D_z)^2);
        deltaL_CD = sqrt((Point_C_x - Point_D_x)^2 + (Point_C_y - Point_D_y)^2 + (Point_C_z - Point_D_z)^2);
        deltaL_AD = sqrt((Point_A_x - Point_D_x)^2 + (Point_A_y - Point_D_y)^2 + (Point_A_z - Point_D_z)^2);
        half_C_ABD = 0.5 * (deltaL_AB + deltaL_BD + deltaL_AD);
        half_C_ACD = 0.5 * (deltaL_AC + deltaL_CD + deltaL_AD);
        y_ABD_max = max([Point_A_y Point_B_y Point_D_y]);
        y_ABD_min = min([Point_A_y Point_B_y Point_D_y]);
        y_ACD_max = max([Point_A_y Point_C_y Point_D_y]);
        y_ACD_min = min([Point_A_y Point_C_y Point_D_y]);
        % 求离散后湿表面积
        S_Triangle_ABD = sqrt(half_C_ABD * (half_C_ABD - deltaL_AB) * (half_C_ABD - deltaL_AD) * (half_C_ABD - deltaL_BD));
        S_Triangle_ACD = sqrt(half_C_ACD * (half_C_ACD - deltaL_AC) * (half_C_ACD - deltaL_AD) * (half_C_ACD - deltaL_CD));
        S_Quadrilateral = S_Triangle_ABD + S_Triangle_ACD;
        S_Discrete_half = S_Discrete_half + S_Quadrilateral; 
    end
end
S_Discrete = 2 * S_Discrete_half