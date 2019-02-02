clc;
clear all;
close all;



%% 1
% 
% M_a = 68.9;                             % Nm
% T_a = 45.7;                             % Nm
% M_m = 57.1;                             % Nm
% T_m = 36.7;                             % Nm
% 
% S_u = 705E6;                            % Pa
% S_y = 565E6;                            % Pa
% 
% S_e = 260E6;                            % Pa
% 
% k_f  = 2.15;
% k_fs = 1.89;
% 
% n = 1.9;
% 
% % DE-Gerber
% A = sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2)
% B = sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)
% d = (8*n*A/(pi*S_e) * (1 + sqrt(1 + (2*B*S_e/(A*S_u))^2)))^(1/3);
% 
% % DE-ASME Elliptic
% d = (16*n/pi * sqrt(4*(k_f*M_a/S_e)^2 + 3*(k_fs*T_a/S_e)^2 ...
%                   + 4*(k_f*M_m/S_y)^2 + 3*(k_fs*T_m/S_y)^2))^(1/3);
% 
% % DE_Soderberg
% d = (16*n/pi * (1/S_e * sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2) ...
%              +  1/S_y * sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)))^(1/3);
% 
% % DE-Goodman
% d = (16*n/pi * (1/S_e * sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2) ...
%              +  1/S_u * sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)))^(1/3);

%% 2
% 
% M_a = 5317;                             % in*lb
% M_c = 6750;                             % in*lb
% T   = 2819;                             % in*lb
% 
% S_u = 63.8;                             % ksi
% S_y = 53.7;                             % ksi
% 
% Se_p = S_u/2;                           % ksi
% if Se_p > 100
%     Se_p = 100;                         % ksi
% end
% 
% % Analyze point C!
% d = 1.3;
% 
% a = 2.7;
% b = -0.265;
% k_a = a*S_u^b;
% 
% k_b = (d/0.3)^-0.107;
% k_c = 1;
% k_d = 1;
% k_e = 1;
% 
% S_e = Se_p * k_a*k_b*k_c*k_d*k_e;       % ksi
% 
% q = 0.6
% k_t = 2.0
% k_ts = 1.7
% k_f = 1 + q*(k_t-1);
% k_fs = 1 + q*(k_ts-1);
% 
% % DE-Gerber
% M_a = M_c;
% T_a = 0;
% M_m = 0;
% T_m = T;
% S_u = S_u * 1000;                       % psi
% S_y = S_y * 1000;                       % psi
% S_e = S_e * 1000;                       % psi
% A = sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2)
% B = sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)
% n = (8*A/(pi*d^3*S_e) * (1 + sqrt(1 + (2*B*S_e/(A*S_u))^2)))^-1

%% 3
% 
% omega = 850;                            % rpm
% n = 1.5;
% 
% S_u = 120;                              % ksi
% S_y = 66;                               % ksi
% 
% Se_p = S_u/2;                           % ksi
% if Se_p > 100
%     Se_p = 100;                         % ksi
% end
% 
% a = 2.7;
% b = -0.265;
% k_a = a*S_u^b;
% 
% k_b = 1;
% k_c = 1;
% k_d = 1;
% k_e = 1;
% 
% S_e = Se_p * k_a*k_b*k_c*k_d*k_e;       % ksi
% 
% f = 0.82;
% a = (f*S_u)^2 / S_e;
% b = -1/3 * log10(f*S_u/S_e);
% N = omega*60*10;                        % cycles (10 hours)
% S_f = a*N^b;
% 
% k_t = 1.65;
% k_ts = 0;
% k_f = k_t;
% k_fs = 0;
% 
% % DE-Goodman
% M_a = 30E3;                             % in*lb
% T_a = 0;
% M_m = 0;
% T_m = 0;
% S_u = S_u * 1000;                       % psi
% S_y = S_y * 1000;                       % psi
% S_f = S_f * 1000;                       % psi
% d = (16*n/pi * (1/S_f * sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2) ...
%              +  1/S_u * sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)))^(1/3)

%% Supplemental 1
A = [-sind(25),      1, 0, 1,    0;
     cosd(25),       0, 1, 0,    1;
     0.15*cosd(25),  0, 0, 0,    0;
     0.75*cosd(25),  0, 0, 0,    1.05;
     -0.75*sind(25), 0, 0, 1.05, 0];
B = [11*sind(20);
     11*cosd(20);
     3.3*cosd(20);
     4.4*cosd(20);
     4.4*sind(20)];
X = linsolve(A, B)