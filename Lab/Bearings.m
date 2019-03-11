clc;
clear all;
close all;



%% Known
L_R = 10^6;

x_0 = 0.02;
theta = 4.459;
b = 1.483;

R_D = 0.999;

%% Carrier Bearing (left)

%% Carrier Bearing (right)

%% Spider Bearing (x2)

%% Middle Drive Bearing (left)

% Known
a = 3;                                  % Ball Bearing
a_f = 1;
F_D = 1;
L_D = 1;

% Analysis
x_D = L_D / L_R;
C_10 = a_f*F_D * (x_D / (x_0 + (theta-x_0) * log(1/R_D)^(1/b)))^(1/a);