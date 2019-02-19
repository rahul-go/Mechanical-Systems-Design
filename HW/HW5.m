clc;
clear all;
close all;



%% 1

% % Given
% d = 1;                                  % in
% d_tol = -0.0015;                        % in
% b = 1.0015;                             % in
% b_tol = 0.003;                          % in
% l_d = 1;
% W = 250;                                % lb
% N = 1100 / 60;                          % rps
% mu = 8E-6;                              % reyn
% 
% % Analysis
% r = d / 2;                              % in
% c = b/2 - r;                            % in
% l = l_d * d;                            % in
% P = W / (2*r*l);                        % psi
% S = (r/c)^2 * mu*N/P
% 
% % From Figure 12-16
% h0_c = 0.59;
% h_0 = h0_c * c                          % in
% 
% % From Figure 12-18
% r_c_f = 5.5;
% f = r_c_f / (r/c);
% 
% H = N * f * W * r;                      % rev/s*in*lb
% H = H * 2*pi / (778*12)                 % Btu/s
% 
% % From Figure 12-19
% Q_rcNl = 3.95;
% Q = Q_rcNl * r*c*N*l;                   % in^3/s
% 
% % From Figure 12-20
% Qs_Q = 0.5;
% Q_s = Qs_Q * Q                          % in^3/s



%% 2

% % Given
% d = 2.580;                              % in
% d_tol = -0.001;                         % in
% b = 2.585;                              % in
% b_tol = 0.004;                          % in
% l = 1.290;                              % in
% N = 600 / 60;                           % rps
% W = 853;                                % lb
% T = 150;                                % degF
% 
% % Analysis
% r = d / 2;                              % in
% c = b/2 - r;                            % in
% P = W / (2*r*l);                        % psi
% l_d = l / d
% 
% 
% 
% % SAE 10
% % From Figure 12-12
% mu = 1.7E-6;                            % reyn
% S = (r/c)^2 * mu*N/P
% 
% % From Figure 12-16
% h0_c = 0.065;
% h_0 = h0_c * c                          % in
% 
% % From Figure 12-21
% p_pmax = 0.175;
% p_max = (p_pmax / P)^-1                 % psi
% 
% 
% 
% % SAE 40
% % From Figure 12-12
% mu = 4.5E-6;                            % reyn
% S = (r/c)^2 * mu*N/P
% 
% % From Figure 12-16
% h0_c = 0.13;
% h_0 = h0_c * c                          % in
% 
% % From Figure 12-21
% p_pmax = 0.23;
% p_max = (p_pmax / P)^-1                 % psi



%% 3

% % Given
% d = 88.00;                              % mm
% d_tol = -0.03;                          % mm
% b = 88.10;                              % mm
% b_tol = 0.06;                           % mm
% l = 44;                                 % mm
% W = 4 * 1000;                           % N
% N = 780 / 60;                           % rps
% T = 60;                                 % degC
% 
% % Analysis
% r = d / 2;                              % mm
% c = b/2 - r;                            % mm
% P = W / (2*r/1000*l/1000);              % Pa
% l_d = l / d
% 
% 
% 
% % SAE 20
% % From Figure 12-13
% mu = 18E-3;                             % Pa*s
% S = (r/c)^2 * mu*N/P
% 
% % From Figure 12-16
% h0_c = 0.29;
% h = h0_c * c                            % in
% 
% % From Figure 12-18
% r_c_f = 5.5;
% f = r_c_f / (r/c);
% 
% H = N * f * W * r;                      % rev/s*N*mm
% H = H * 2*pi/1000                       % W
% 
% % From Figure 12-21
% p_pmax = 0.32;
% p_max = (p_pmax / P)^-1
% 
% 
% 
% % SAE 40
% % From Figure 12-13
% mu = 37E-3;                             % Pa*s
% S = (r/c)^2 * mu*N/P
% 
% % From Figure 12-16
% h0_c = 0.44;
% h = h0_c * c                            % in
% 
% % From Figure 12-18
% r_c_f = 8.75;
% f = r_c_f / (r/c);
% 
% H = N * f * W * r;                      % rev/s*N*mm
% H = H * 2*pi/1000                       % W
% 
% % From Figure 12-21
% p_pmax = 0.38;
% p_max = (p_pmax / P)^-1



%% 4

% Given
l = 1;                                  % in
d = 1;                                  % in
T = 70;                                 % degF
w = 0.005;                              % in
W = 572;                                % lb
N = 222;                                % rpm

% Analysis
P = 4 / pi * W / (d * l)                % psi
V = N * pi*d/12                         % ft/min

% From Table 12-8
K = 0.6E-10;
PV = 46700;

% From Table 12-9
f_s = 0.03;

% From Table 12-10
% f_1 = 1.3 + (V - 33)/(100 - 33)*(1.8-1.3);
f_1 = 1.8

% From Table 12-11
f_2 = 1.0;

if P*V < PV
    PV = P*V;
end

% w = f_1 * f_2 * K * P * V * t           % in
t = w / (f_1 * f_2 * K * PV)            % hr
N * 60 * t