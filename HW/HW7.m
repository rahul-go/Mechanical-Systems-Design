clc;
clear all;
close all;



%% 1

% % Given
% d = 7;                                  % mm
% C = 10;
% L = 81;                                 % mm
% F = 44;                                 % N
% delta = 14;                             % mm
% 
% % Analysis
% k = F / delta                           % N/mm
% 
% % From Table 10-5
% G = 77.2E3;                             % kPa
% 
% D = d * C;                              % in
% N = (d)^4*G / (8*(D)^3*k)               % coils
% N = (d)^4*G / (8*(D)^3*k) / (1 + 1/(2*C^2))



%% 2

% % Given
% d = 0.072;                              % in
% OD = 0.896;                             % in
% N_t = 8;
% 
% % Analysis
% D = OD - d;
% C = D / d;
% n_s = 1.2;
% 
% % From Table 10-4
% A = 140;
% m = 0.190;
% S_ut = A / d^m;                         % ksi
% 
% % From Table 10-6
% S_sy = 0.45 * S_ut;                     % ksi
% 
% tau = S_sy * 1000;                             % ksi
% K_b = (4*C + 2) / (4*C - 3);
% F = tau/n_s * pi*d^3 / (K_b * 8*D);     % lb
% 
% % From Table 10-5
% G = 11.5E6;                             % psi
% 
% N_a = N_t - 1;
% k = d^4*G / (8*D^3*N_a);
% k = d^4*G / (8*D^3*N_a) / (1 + 1/(2*C^2));
% 
% y = F / k;
% 
% L_s = d * N_t;
% L_0 = L_s + y
% 
% p = L_0 / N_t
% 
% alpha = S_sy / n_s;
% L0_cr = 2.63 * D / 0.5
% L_0 < L0_cr



%% 3

% % Given
% d = 0.051;                              % in
% OD = 0.255;                             % in
% L_0 = 0.68;                             % in
% N_t = 11.9;
% 
% % Analysis
% D = OD - d;                             % in
% C = D / d;
% 
% % From Table 10-4
% A = 169;
% m = 0.146;
% S_ut = A / d^m * 10^3;                  % psi
% 
% % From Table 10-6
% S_sy = 0.35 * S_ut;                     % psi
% 
% % From Table 10-5
% G = 10E6;                               % psi
% 
% N_a = N_t - 2;
% k = d^4*G / (8*D^3*N_a);
% L_s = d * N_t;
% y = L_0 - L_s;
% F = k * y;
% 
% K_b = (4*C + 2) / (4*C - 3);
% tau = K_b * 8*F*D / (pi*d^3);
% 
% n_s = S_sy / tau



%% 4

% Given
t = 37.5;                               % mm
F = 55;                                 % N
screw_OD = 10;                          % mm
d_c = 1.25;                             % mm
L0_max = 48;                            % mm
Ls_max = 31.5;                          % mm
ns_min = 1.2;
d = 2;                                  % mm



% Pre-Analysis
ID_min = screw_OD + d_c;                % mm

% % Iteration 1
d = d;                                  % mm
D = ID_min + d;                         % mm
L_s = Ls_max;                           % mm
L_0 = L0_max;
OD = ID_min + 2*d;                      % mm



% Analysis
C = D / d;

% From Table 10-4
A = 1783;
m = 0.190;
S_ut = A / d^m;                         % MPa

% From Table 10-5
G = 79.3E3;                             % MPa

% From Table 10-6
S_sy = 0.45 * S_ut;                     % MPa

y = L_0 - t;
k = F / y;
N_a = d^4*G / (8*D^3*k);
N_t = N_a + 2;

K_b = (4*C + 2) / (4*C - 3);
y_s = L_0 - L_s;                        % mm
F_s = k * y_s;                          % N
tau_s = K_b * 8*F_s*D / (pi*d^3);       % MPa

% FOS
S_sy / tau_s