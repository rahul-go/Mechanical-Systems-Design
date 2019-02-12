clc;
clear all;
close all;



%% 1

% % Given
% a = 10/3;
% F_D = 16;                               % kN
% hours = 8000;                           % hr
% speed = 920;                            % rpm
% R = 0.95;
% L_R = 10^6;                             % rev
% x_0 = 0.02;
% theta = 4.459;
% b = 1.483;
% 
% % Analysis
% L_D = speed * hours * 60;
% x_D = L_D / L_R;
% C_10 = F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)
% 
% % Analysis
% L_D = speed * hours * 60;
% C_10 = F_D * (L_D/L_R)^(1/a)



%% 2

% % Given
% a = 3;
% L_R = 10^6;                             % rev
% F_D = 4.1;                              % kN
% L_D = 10^9;
% R = 0.90;
% 
% % Analysis
% C_10 = F_D * (L_D/L_R)^(1/a)



%% 3

% % Given
% a = 3;
% L_R = 10^6;                             % rev
% a_f = 1;
% x_0 = 0.02;
% theta = 4.459;
% b = 1.483;
% F_D = 644;                              % lb
% years = 5;                              % years
% hrperwk = 40;                           % hr/week
% speed = 393;                            % rpm
% R = 0.95;
% 
% % Analysis
% L_D = years * 52 * hrperwk * 60 * speed;
% x_D = L_D / L_R;
% C_10 = a_f * F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)



%% 4

% % Given
% V = 1;
% bore = 30;                              % mm
% F_a = 2;                                % kN
% F_r = 4;                                % kN
% speed = 402;                            % rpm
% 
% % Analysis
% C_0 = 10;
% F_a / C_0
% e = (0.34*8 + 0.38*3) / 11;
% if (F_a/(V*F_r) <= e)
%     % 1
%     X = 1;
%     Y = 0;
% else
%     % 2
%     X = 0.56;
%     Y = (1.31*8 + 1.15*3) / 11;
% end
% F_e = X*V*F_r + Y*F_a



%% 5

% % Given
% a = 3;
% L_R = 10^6;                             % rev
% a_f = 1;
% x_0 = 0.02;
% theta = 4.459;
% b = 1.483;
% F_r = 9;                                % kN
% F_a = 5;                                % kN
% hours = 12000;                          % hr
% speed = 250;                            % rpm
% V = 1;
% R = 0.95;

% Analysis

% % C_0
% % Start from middle line
% X = 0.56;
% Y = 1.63;
% F_e = X*V*F_r + Y*F_a
% 
% % C_10
% % Find associated bearing
% F_D = F_e;
% L_D = speed * hours * 60;
% x_D = L_D / L_R;
% C_10 = a_f * F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)
% % 90 mm bearing
% 
% % C_0
% % Iterate
% C_0 = 62.0;
% F_a / C_0
% e = 0.27;
% F_a / (V * F_r)
% % Checks out.



% % TRY DIFFERENT STARTING POINT
% % C_0
% % Start from last line
% X = 0.56;
% Y = 1.00;
% F_e = X*V*F_r + Y*F_a;
% 
% % C_10
% % Find associated bearing
% F_D = F_e;
% L_D = speed * hours * 60;
% x_D = L_D / L_R;
% C_10 = a_f * F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a);
% % 80 mm bearing
% 
% % C_0
% % Iterate 1
% C_0 = 45.0;
% F_a / C_0;
% e = 0.30;
% F_a / (V * F_r);
% X = 0.56;
% Y = 1.45;
% F_e = X*V*F_r + Y*F_a;
% 
% % C_10
% % Iterate 1
% % Find associated bearing
% F_D = F_e;
% L_D = speed * hours * 60;
% x_D = L_D / L_R;
% C_10 = a_f * F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a);
% % 85 mm bearing
% 
% % C_0
% % Iterate 2
% C_0 = 53.0;
% F_a / C_0;
% e = 0.29;
% F_a / (V * F_r);
% X = 0.56;
% Y = 1.50;
% F_e = X*V*F_r + Y*F_a;
% 
% % C_10
% % Iterate 2
% % Find associated bearing
% F_D = F_e;
% L_D = speed * hours * 60;
% x_D = L_D / L_R;
% C_10 = a_f * F_D * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)
% % 85 mm bearing
% % Checks out.



%% 6

% Given
a = 10/3;
speed = 375;                            % rpm
hrperday = 8;                           % hr/day
dayperwk = 5;                           % day/week
years = 5;                              % years
R = 0.90;
F_ru = 14;                              % kN
F_rl = 33;                              % kN
F_a = 6;                                % kN
K = 1.5;
a_f = 1.2;
L_R = 90*10^6;

% Analysis
x_0 = 0;
theta = 4.48;
b = 3/2;

L_D = speed * 60 * hrperday * dayperwk * 52 * years;
x_D = L_D / L_R;

% Lower: Bearing A | Upper: Bearing B
F_iu = 0.47 * F_ru / K;
F_il = 0.47 * F_rl / K;
F_il;
F_iu + F_a;
% Fia <= (FiB + Fae) --> First equations

% Upper Bearing
F_eu = F_ru;
F_Du = F_eu;
C_10u = a_f * F_Du * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)

% % Lower Bearing
F_el = 0.4*F_rl + K*(F_iu+F_a);
% Use radial load, F_rl for calculations
F_Dl = F_rl;
C_10l = a_f * F_Dl * (x_D / (x_0 + (theta-x_0)*log(1/R)^(1/b)))^(1/a)