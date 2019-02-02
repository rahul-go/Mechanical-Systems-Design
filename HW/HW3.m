clc;
clear all;
close all;



% %% 1
% 
% % Given
% phi = 20;                               % deg
% N_p = 20;                               % teeth
% H = 120;                                % watts
% N_g = 36;                               % teeth
% m = 3.0;                                % mm
% omega_p = 90 * 2*pi/60;                 % rad/s
% b = 18;                                 % mm
% H_B = 200;
% N_L = 10^8;                             % cycles
% R = 0.95;
% 
% % Analysis
% d_p = N_p * m;                          % mm
% T_p = H / omega_p;                      % N*m
% d_g = d_p * N_g/N_p;                    % mm
% W_t = T_p / (d_p/2/1000);               % N
% 
% %% Pinion Stresses
% Y_p = 0.322;
% % Overload Factor
% K_o = 1;
% % Dynamic Factor
% V = omega_p * d_p/2/1000;               % m/s
% B = 0.25 * (12-6)^(2/3);
% A = 50 + 56*(1-B);
% K_v = ((A + sqrt(200*V)) / A)^(B);
% % Size Factor
% K_s = 0.8433 * (m*b*sqrt(Y_p))^0.0535;  % Connect
% % K_s = 1.192 * (m_p*b*sqrt(Y_p))^0.0535; % Shigley's
% % Load-Distribution Factor
% C_mc = 1;
% C_pf = b / (10*d_p) - 0.025;
% if b / (10*d_p) < 0.05
%     C_pf = 0.025;
% end
% C_pm = 1;
% A = 0.247;
% B = 0.0167;
% C = -0.765*10^-4;
% C_ma = A + B*(b/25.4) + C*(b/25.4)^2;
% C_e = 1;
% K_H = 1 + C_mc * (C_pf*C_pm + C_ma*C_e);
% % Rim-Thickness Factor
% K_B = 1;
% % Bending-Strength Geometry Factor
% Y_J = 0.325;
% % Bending Stress
% sigma_p = W_t*K_o*K_v*K_s/(b*m)*K_H*K_B/Y_J % MPa
% % Elastic Coefficient
% Z_E = 191;                              % MPa^(1/2)
% % Surface Condition Factor
% Z_R = 1;
% % Surface-Strength Geometry Factor
% m_G = N_g / N_p;
% m_N = 1;
% Z_I = cosd(20)*sind(20)/(2*m_N) * m_G/(m_G+1);
% % Contact Stress
% sigma_pc = Z_E * sqrt(W_t*K_o*K_v*K_s*K_H/(d_p*b)*Z_R/Z_I) % MPa
% 
% %% Gear Stresses
% Y_g = (0.371 + 0.384) / 2;
% % Overload Factor
% K_o = 1;
% % Dynamic Factor
% V = V;                                  % m/s
% B = 0.25 * (12-6)^(2/3);
% A = 50 + 56*(1-B);
% K_v = ((A + sqrt(200*V)) / A)^(B);
% % Size Factor
% K_s = 0.8433 * (m*b*sqrt(Y_g))^0.0535;  % Connect
% % K_s = 1.192 * (m_p*b*sqrt(Y_g))^0.0535; % Shigley's
% % Load-Distribution Factor
% C_mc = 1;
% C_pf = b / (10*d_g) - 0.025;
% if b / (10*d_g) < 0.05
%     C_pf = 0.025;
% end
% C_pm = 1;
% A = 0.247;
% B = 0.0167;
% C = -0.765*10^-4;
% C_ma = A + B*(b/25.4) + C*(b/25.4)^2;
% C_e = 1;
% K_H = 1 + C_mc * (C_pf*C_pm + C_ma*C_e);
% % Rim-Thickness Factor
% K_B = 1;
% % Bending-Strength Geometry Factor
% Y_J = 0.375;
% % Bending Stress
% sigma_g = W_t*K_o*K_v*K_s/(b*m)*K_H*K_B/Y_J % MPa
% % Elastic Coefficient
% Z_E = 191;                              % MPa^(1/2)
% % Surface Condition Factor
% Z_R = 1;
% % Surface-Strength Geometry Factor
% m_G = m_G;
% m_N = 1;
% Z_I = cosd(20)*sind(20)/(2*m_N) * m_G/(m_G+1);
% % Contact Stress
% sigma_gc = Z_E * sqrt(W_t*K_o*K_v*K_s*K_H/(d_p*b)*Z_R/Z_I) % MPa
% 
% %% Allowable Stresses
% % Bending Stress
% S_t = 0.533 * H_B + 88.3;               % MPa
% % Stress-Cycle Factor
% Y_N = 1.3558*N_L^-0.0178;
% % Y_N = (1.3558*N_L^-0.0178 + 1.6831*R^-0.0323) / 2;
% % Temperature Factor
% Y_theta = 1;
% % Reliability Factor
% Y_Z = 0.658 - 0.0759*log(1 - R);
% % Bending Stress
% sigma_all = S_t/1 * Y_N/(Y_theta*Y_Z);  % MPa
% % Contact Stress
% S_c = 2.22 * H_B + 200;                 % MPa
% % Stress-Cycle Factor
% Z_N_p = 1.4488*N_L^-0.023;
% Z_N_g = 1.4488*(N_L/(N_g/N_p))^-0.023;
% % Hardness Ratio Factor
% Z_W_p = 1;
% A_prime = 0;
% Z_W_g = 1.0 + A_prime*(m_G - 1.0);
% % Temperature Factor
% Y_theta = Y_theta;
% % Reliability Factor
% Y_Z = Y_Z;
% % Contact Stress
% sigma_all_cp = S_c * Z_N_p*Z_W_p/(Y_theta*Y_Z); % MPa
% sigma_all_cg = S_c * Z_N_g*Z_W_g/(Y_theta*Y_Z); % MPa
% 
% %% FOS
% FOS_p = sigma_all / sigma_p
% FOS_cp = sigma_all_cp / sigma_pc
% FOS_g = sigma_all / sigma_g
% FOS_cg = sigma_all_cg / sigma_gc



% %% 2
% 
% % Given
% N_p = 20;                               % teeth
% P_d = 6;                                % teeth/in
% Q_v = 6;
% H_B = 300;
% N_g = 60;                               % teeth
% N_L = 10^9;                             % cycles
% R = 0.999;
% angle = 90;                             % deg
% omega_p = 750 * 2*pi/60;                % rad/s
% F = 1.20;                               % in
% phi = 20;                               % deg
% K_0 = 1;
% S_F = 1;
% S_H = 1;
% 
% % Analysis
% d_p = 20 / P_d;
% d_g = 60 / P_d;
% omega_g = omega_p / (N_g/N_p);
% 
% %% Allowable Pinion Bending Stress
% s_at = 44 * H_B + 2100;
% K_L = 1.683*N_L^-0.0323;
% K_T = 1;
% K_R = 0.50 - 0.25*log10(1-R);
% s_wt = s_at*K_L / (S_F*K_T*K_R);
% 
% % Pinion Bending Stress
% s_t = s_wt;
% % Overload Factor
% K_o = 1;
% % Dynamic Factor
% v_t = omega_p * d_p/2/12 * 60;          % ft/min
% B = 0.25 * (12 - Q_v)^(2/3);
% A = 50 + 56 * (1 - B);
% K_v = ((A + sqrt(v_t)) / A)^B;
% % Size Factor
% K_s = 0.4867 + 0.2132 / P_d;
% % Load-Distribution Factor
% K_mb = 1.1;
% K_m = K_mb + 0.0036*F^2;
% % Lengthwise Curvature Factor
% K_x = 1;
% % Bending Strength Geometry Factor
% J = 0.248;
% W_t = s_t * (1/F*P_d*K_o*K_v*K_s*K_m/(K_x*J))^-1;
% 
% % Rated Power
% H_p = W_t * v_t / 33000                 % hp
% 
% %% Allowable Gear Bending Stress
% s_at = 44 * H_B + 2100;
% K_L = 1.683*(N_L / (N_g/N_p))^-0.0323;
% K_T = 1;
% K_R = 0.50 - 0.25*log10(1-R);
% s_wt = s_at*K_L / (S_F*K_T*K_R);
% 
% % Gear Bending Stress
% s_t = s_wt;
% % Overload Factor
% K_o = 1;
% % Dynamic Factor
% v_t = v_t;
% B = 0.25 * (12 - Q_v)^(2/3);
% A = 50 + 56 * (1 - B);
% K_v = ((A + sqrt(v_t)) / A)^B;
% % Size Factor
% K_s = 0.4867 + 0.2132 / P_d;
% % Load-Distribution Factor
% K_mb = 1.1;
% K_m = K_mb + 0.0036*F^2;
% % Lengthwise Curvature Factor
% K_x = 1;
% % Bending Strength Geometry Factor
% J = 0.202;
% W_t = s_t * (1/F*P_d*K_o*K_v*K_s*K_m/(K_x*J))^-1;
% 
% % Rated Power
% H_g = W_t * v_t / 33000                 % hp