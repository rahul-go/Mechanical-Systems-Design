clc;
clear all;
close all;



%% 2

% % Given
% pitch = 5;                              % teeth/in
% teeth_2 = 16;
% teeth_3 = 40;
% H = 33 * 550*12;                        % in*lb/s
% omega = 1800 * 2*pi/60;                 % rad/s
% 
% % Analysis
% d_2 = teeth_2 / pitch;                  % in
% r_2 = d_2 / 2;                          % in
% d_3 = teeth_3 / pitch;                  % in
% r_3 = d_3 / 2;                          % in
% 
% T = H/omega;                            % in*lb
% F_2z = T/r_2;
% F_2y = -F_2z*tand(20);
% A_y = -F_2y / 2;
% A_z = -F_2z / 2;
% R_A = hypot(A_y, A_z)
% B_y = -F_2y / 2;
% B_z = -F_2z / 2;
% R_B = hypot(B_y, B_z)
% 
% F_3z = -F_2z;
% F_3y = -F_2y;
% C_y = -F_3y / 2;
% C_z = -F_3z / 2;
% R_C = hypot(C_y, C_z)
% D_y = -F_3y / 2;
% D_z = -F_3z / 2;
% R_D = hypot(D_y, D_z)



%% 3

% % Given
% phi_n = 20;                             % deg
% psi = 30;                               % deg
% W_3t = 359;                             % lb
% pitch = 7;                              % teeth/in
% teeth_3 = 54;
% teeth_4 = 14;
% 
% % Analysis
% r_3 = teeth_3 / pitch;                  % in
% r_4 = teeth_4 / pitch;                  % in
% 
% W_3r = W_3t * tand(phi_n) / cosd(psi)   % lb
% w_3a = W_3t * tand(psi)                 % lb
% 
% W_4t = W_3t * (r_3/r_4)                 % lb
% W_4r = W_4t * tand(phi_n) / cosd(psi)   % lb
% W_4a = W_4t * tand(psi)                 % lb

%% 4

% Given
H = 1700;                               % watts
omega = 500 * 2*pi/60;                  % rad/s
pitch = 25 / 1000;                      % m/tooth
phi_n = 14.5;                           % deg
d_w = 100 / 1000;                       % m
w_w = 100 / 1000;                       % m
w_g = 50 / 1000;                        % m
AB = 200 / 1000;                        % m

% Analysis
T = H / omega                           % Nm

V_w = omega * d_w/2;                    % m/s
W_x = T / (d_w/2)                       % N
L = pitch;                              % m
gamma = atand(L/(pi*w_w));              % deg
V_s = V_w/cosd(gamma) * 196.85;         % ft/min
f = 0.048;

W = W_x / (cosd(phi_n) * sind(gamma) + f*cosd(gamma));  % N
W_y = W * sind(phi_n)                                   % N
W_z = W * (cosd(phi_n) * cosd(gamma) - f*sind(gamma))   % N

A_x = -W_x / 2                          % N
B_x = -W_x / 2                          % N
A_y = -W_y/2 - W_z/4                    % N
B_y = -W_y/2 + W_z/4                    % N
A_z = -W_z                              % N
B_z = 0                                 % N