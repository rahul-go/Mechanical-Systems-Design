clc;
clear all;
close all;



%% 1

% % Given
% a = 10/3;
% F_D = 16;                               % kN
% hours = 9000;                           % hr
% speed = 920;                            % rpm
% R = 0.95;
% L_R = 10^6;
% x_0 = 0.02;
% theta = 4.459;
% b = 1.483;
% 
% % Analysis
% L_D = speed * hours * 60;
% C_10 = F_D * (L_D/L_R)^(1/a)



%% 2

% % Given
% a = 3;
% L_R = 10^6;
% F_D = 4.1;                              % kN
% L_D = 10^9;
% R = 0.90;
% 
% % Analysis
% C_10 = F_D * (L_D/L_R)^(1/a)



%% 3

% Given
a = 3;
L_R = 10^6;
a_f = 1;
x_0 = 0.02;
theta = 4.459;
b = 1.483;
F_D = 663;                              % lb
years = 5;                              % years
hrperwk = 40;                           % hr/week
speed = 393;                            % rpm
R = 0.95;

% Analysis
L_D = years * 52 * hrperwk * 60 * speed;
C_10 = F_D * (L_D/L_R)^(1/a)