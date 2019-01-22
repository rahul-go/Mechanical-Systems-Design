clc;
clear all;
close all;

%% 2

% Given
pitch = 5;                              % teeth/in
teeth_2 = 16;
teeth_3 = 40;
H = 33 * 550*12;                        % in*lb/s
omega = 1800 * 2*pi/60;                 % rad/s

% Analysis
d_2 = teeth_2 / pitch;                  % in
r_2 = d_2 / 2;                          % in
d_3 = teeth_3 / pitch;                  % in
r_3 = d_3 / 2;                          % in

T = H/omega;                            % in*lb
F_2z = T/r_2;
F_2y = -F_2z*tand(20);
A_y = -F_2y / 2;
A_z = -F_2z / 2;
R_A = hypot(A_y, A_z)
B_y = -F_2y / 2;
B_z = -F_2z / 2;
R_B = hypot(B_y, B_z)

F_3z = -F_2z;
F_3y = -F_2y;
C_y = -F_3y / 2;
C_z = -F_3z / 2;
R_C = hypot(C_y, C_z)
D_y = -F_3y / 2;
D_z = -F_3z / 2;
R_D = hypot(D_y, D_z)