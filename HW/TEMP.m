
M_a = 5320;                             % in*lb
M_c = 6750;                             % in*lb
T   = 2819;                             % in*lb

S_u = 64;                             % ksi
S_y = 55;                             % ksi

Se_p = S_u/2;                           % ksi
if Se_p > 100
    Se_p = 100;                         % ksi
end

% Analyze point C!
d = 1.3;

a = 2.7;
b = -0.265;
k_a = a*S_u^b;

k_b = (d/0.3)^-0.107;
k_c = 1;
k_d = 1;
k_e = 1;

S_e = Se_p * k_a*k_b*k_c*k_d*k_e;       % ksi

q = 0.6
k_t = 2.0
k_ts = 1.7
k_f = 1 + q*(k_t-1);
k_fs = 1 + q*(k_ts-1);

% DE-Gerber
M_a = M_c;
T_a = 0;
M_m = 0;
T_m = T;
S_u = S_u * 1000;                       % psi
S_y = S_y * 1000;                       % psi
S_e = S_e * 1000;                       % psi
A = sqrt(4*(k_f*M_a)^2 + 3*(k_fs*T_a)^2)
B = sqrt(4*(k_f*M_m)^2 + 3*(k_fs*T_m)^2)
n = (8*A/(pi*d^3*S_e) * (1 + sqrt(1 + (2*B*S_e/(A*S_u))^2)))^-1
