clc;
clear all;
close all;
% Review (not HW)



% %% 1
% 
% % Given
% d = 340;                                % mm
% R = 125;                                % mm
% b = 40;                                 % mm
% F = 2.3;                                % kN
% f = 0.26;
% 
% % Analysis
% r = d / 2;                              % mm
% a = R;                                  % mm
% c = 2*R*sind(60);                       % mm
% theta_1 = 0;                            % deg
% theta_2 = 120;                          % deg
% 
% % Guess P_a (kPa) [self-energizing]
% P_a = 705.1;
% P_a = P_a * 1E-3;                       % MPa
% Pa_1 = P_a;                             % MPa
% 
% M_f = f*P_a*b*r/sin(pi/2) * (-r*(cosd(theta_2)-cosd(theta_1)) - a/2*((sind(theta_2))^2-(sind(theta_1))^2));
% M_N = P_a*b*r*a/sin(pi/2) * ((theta_2-theta_1)/2*pi/180 - sind(2*(theta_2-theta_1))/4);
% F = (M_N - M_f) / c * 1E-3              % kN
% 
% % Guess P_a (kPa) [self-deenergizing]
% P_a = 345.3;
% P_a = P_a * 1E-3;                       % MPa
% Pa_2 = P_a;
% 
% M_f = f*P_a*b*r/sin(pi/2) * (-r*(cosd(theta_2)-cosd(theta_1)) - a/2*((sind(theta_2))^2-(sind(theta_1))^2));
% M_N = P_a*b*r*a/sin(pi/2) * ((theta_2-theta_1)/2*pi/180 - sind(2*(theta_2-theta_1))/4);
% F = (M_N + M_f) / c * 1E-3              % kN
% 
% % Torque [self-energizing]
% P_a = Pa_1;                             % MPa
% T = f*P_a*b*r^2/sin(pi/2) * (cosd(theta_1) - cosd(theta_2));
% T = T * 1E-3                            % N*m
% T_1 = T;
% % Torque [self-deenergizing]
% P_a = Pa_2;                             % MPa
% T = f*P_a*b*r^2/sin(pi/2) * (cosd(theta_1) - cosd(theta_2));
% T = T * 1E-3                            % N*m
% T_2 = T;
% % Torque
% T = T_1 + T_2

% Unfinished



%% 2

% Given
f = 0.31;
b = 2;                                  % in
Pa_max = 140;                           % psi

% Analysis
P_a = Pa_max;                           % psi
r = 10;                                 % in
a = sqrt(3^2 + 12^2);                   % in
c_L = 12 + 12 + 4;                      % in
c_R = (24 - 2*tand(14)) * cosd(14);     % in
theta_1 = 6;                            % deg
theta_2 = 136;                          % deg

M_f = f*P_a*b*r/sin(pi/2) * (-r*(cosd(theta_2)-cosd(theta_1)) - a/2*((sind(theta_2))^2-(sind(theta_1))^2));
M_N = P_a*b*r*a/sin(pi/2) * ((theta_2-theta_1)/2*pi/180 - sind(2*(theta_2-theta_1))/4);

% Left (self-energizing)
c = c_L;                                % in
F = (M_N - M_f) / c

% % Right (self-deenergizing)
c = c_R;                                % in
F = (M_N + M_f) / c

% Wrong