clc;
clear all;
close all;



% %% 1
% 
% % Given
% d = 24;                                 % mm
% l = 5;                                  % mm
% F = 3;                                  % kN
% f_c = 0.04;
% f_t = 0.07;
% d_c = 44;                               % mm
% 
% 
% 
% % Analysis
% d_m = (d + (d-l)) / 2;
% T_r = F*d_m/2 * (l + pi*f_t*d_m)/(pi*d_m - f_t*l);
% T_c = F * f_c * d_c / 2;
% T_r + T_c
% 
% T_l = F*d_m/2 * (pi*f_t*d_m - l)/(pi*d_m + f_t*l);
% T_l + T_c
% 
% e = F * l / (2 * pi * (T_r+T_c))



%% 2

% % Given
% N = 6;
% k_b = 4;                                % Mlb/in
% k_m = 16;                               % Mlb/in
% P_total = 94E3;                         % lb
% d = 1/2;                                % in
% p = 13;                                 % teeth/in
% preload = 0.75;                         % ...of proof load
% 
% % Analysis
% 
% % From Table 8-2
% A_t = 0.1419                            % in^2
% % From Table 8-9
% S_p = 120E3;                            % ksi
% S_t = 150E3;                            % ksi
% S_y = 130E3;                            % ksi
% % -----
% F_p = A_t * S_p;
% F_t = A_t * S_t;
% F_y = A_t * S_y;
% 
% C = k_b / (k_b + k_m);
% P = P_total/N;                          % lb
% F_i = preload * F_p;
% P_b = C * P + F_i;
% 
% %%% Learn this! Interesting difference between the two.
% %%% Yielding FOS vs. Overloading FOS AKA Load Factor
% n_p = F_p / P_b
% n_L = (F_p - F_i) / (C * P)



%% 3

% % Given
% N = 10;
% d_s = 151;                              % mm
% A = 100;                                % mm
% B = 200;                                % mm
% C = 300;                                % mm
% D = 20;                                 % mm
% E = 25;                                 % mm
% p = 7;                                  % MPa
% d_b = 12;                               % mm
% 
% % Analysis
% preload = 0.75;                         % ...of proof load
% 
% % From Table 8-1
% A_t = 84.3;                             % mm
% d_w = 18;                               % mm
% 
% P = p * pi*(d_s/2)^2 / N;
% 
% l = D + E;                              % mm
% % From Table A-31
% H = 10.8;                               % mm
% L_t = 2 * d_b + 6;                      % mm
% l + H
% % From Table A-17
% L = 60;                                 % mm
% l_d = L - L_t;                          % mm
% l_t = l - l_d;                          % mm
% 
% A_d = pi * (d_b/2)^2;                   % mm
% 
% % From MAGIC
% E_0 = 207;                              % GPa
% k_b = (A_d*A_t*E_0) / (A_d*l_t + A_t*l_d);
% 
% % Steel Cylinder Head
% t = D;                                  % mm
% D_ = d_w;
% % From MAGIC
% E_ = 207;
% k_1 = 1/sqrt(3) * pi * E_ * d_b / ...   % MN/m
%     log((1.155*t+D_-d_b)*(D_+d_b) / ((1.155*t+D_+d_b)*(D_-d_b)));
% 
% % Upper Frustum
% midpoint = (D + E) / 2;                 % mm
% t = midpoint - D;                       % mm
% D_ = d_w + 2 * D * tand(30);            % mm
% % From MAGIC
% E_ = 100;                               % GPa
% k_2 = 1/sqrt(3) * pi * E_ * d_b / ...   % MN/m
%     log((1.155*t+D_-d_b)*(D_+d_b) / ((1.155*t+D_+d_b)*(D_-d_b)));
% 
% % Lower Frustum
% t = midpoint;                           % mm
% D_ = d_w;                               % mm
% % From MAGIC
% E_ = 100;                               % GPa
% k_3 = 1/sqrt(3) * pi * E_ * d_b / ...   % MN/m
%     log((1.155*t+D_-d_b)*(D_+d_b) / ((1.155*t+D_+d_b)*(D_-d_b)));
% 
% k_m = 1 / (1/k_1 + 1/k_2 + 1/k_3);      % MN/m
% 
% C_ = k_b / (k_b + k_m);
% 
% % From Table 8-11
% S_p = 600;                              % MPa
% F_p = A_t * S_p;
% F_i = preload * F_p;
% 
% n = (S_p*10^-3 * A_t - F_i/10^3) / (C_ * P/10^3)

%% 4

% Given
F_s = 3977;                             % lb
d = 0.375;                              % in
t = 0.25;                               % in
N = 2;

% Analysis
S_yc = 32E3;                            % psi
S_y = 92E3;                             % psi
S_sy = 1/sqrt(3) * S_y;                 % psi

% Shear of bolts
A_s = pi * (d/2)^2 * N;                 % in^2
tau = F_s / A_s;                        % ksi
n = S_sy / tau

% Bearing on bolts
A_b = d * t * N;                        % in^2
sigma_b = F_s / A_b;                    % ksi
n = S_y / sigma_b

% Bearing on members
n = S_yc / sigma_b

% Tension of members
% ???
A_t = (2.375 - 0.75)*t;                 % in^2
sigma_t = F_s / A_t;                    % ksi
n = S_yc / sigma_t



%% 5

% Given
l = 152;                                % mm
w = 76;                                 % mm
t = 12;                                 % mm
d = 12;                                 % mm
p = 1.75;                               % mm
n_d = 2.6;

% Analysis
% TODO