%% Four Tank System_Phase 2 (Non-linear)
clc;
clear;
close all;

%% Estimated parameter values of the real plant
% nominal levels in cm
h1_eq = 11.4;
h2_eq = 11.6;
h3_eq = 5.3;
h4_eq = 4;

% nominal pump settings
v1_eq = 0.5;
v2_eq = 0.5;
a = [2.10 2.14 2.2 2.3];  % area of the drain in cm^2
Area = [730 730 730 730]; % area of the tanks in cm^2
y = [0.3 0.35];           % y(1) = Ratio of flow in tank1 to flow in tank4
                          % y(2) = Ratio of flow in tank2 to flow in tank3                
k = [7.45 7.30];          % pump proportionality constants in cm^2/s
g = 981;                  % gravitational acceleration in cm/s^2
% taw =[2 2.1];             % pump response time constants in s

%% Non-linear differential equations
syms h1 h2 h3 h4 h1_d h2_d h3_d h4_d v1 v2 
h1d = -a(1)/Area(1)*sqrt(2*g*h1)+a(3)/Area(1)*sqrt(2*g*h3)+y(1)*k(1)/Area(1)*v1;
h2d = -a(2)/Area(2)*sqrt(2*g*h2)+a(4)/Area(2)*sqrt(2*g*h4)+y(2)*k(2)/Area(2)*v2;
h3d = -a(3)/Area(3)*sqrt(2*g*h3)+(1 - y(2))*(k(2)/Area(3))*v2;
h4d = -a(4)/Area(4)*sqrt(2*g*h4)+(1 - y(1))*(k(1)/Area(4))*v1;
% v1d = -v1/taw(1)+u1/taw(1);
% v2d = -v2/taw(2)+u2/taw(2);

%% Linearization and Steady-State model
F = [h1d; h2d; h3d; h4d];

E1 = subs(h1d, [h1,h2], [h1_eq, h2_eq]);
E2 = subs(h2d, [h1,h2], [h1_eq, h2_eq]);
E3 = subs(h3d, [h1,h2], [h1_eq, h2_eq]);
E4 = subs(h4d, [h1,h2], [h1_eq, h2_eq]);

eq_values = zeros(4,1);
tmp = vpasolve([E1==0; E2==0; E3==0; E4==0], [h3 h4 v1 v2]);
total_ans = fieldnames(tmp);
for i=1:4
    eq_values(i,1) = real(double(getfield(tmp, total_ans{i})));
end

Amat = jacobian(F,[h1,h2,h3,h4]);
Amat = subs(Amat,[h1,h2,h3,h4,v1,v2],[h1_eq,h2_eq,h3_eq,h4_eq,v1_eq,v2_eq]);
A = double(Amat);

Bmat = jacobian(F,[v1,v2]);
Bmat = subs(Bmat,[h1,h2,h3,h4,v1,v2],[h1_eq,h2_eq,h3_eq,h4_eq,v1_eq,v2_eq]);
B = double(Bmat);

% outputs are h1 & h2 (level of tank 1 and tank 2)
C=[1 0 0 0;
   0 1 0 0];
 
D = [0 0 ;
    0 0];

non_linear_sys = ss(A,B,C,D);

%% Non-linear static compensator
P_non_stat = [-0.2, -0.15, -0.04, -0.08];
K_non_stat = place(A, B, P_non_stat);
K_non_stat = vpa(K_non_stat, 3)
G_cl0_non = vpa(inv(C*inv(-(A-B*K_non_stat))*B), 3)


%% Dynamic Compensator
% check for the necessary condition
if rank([B, A; zeros(2,2), -C]) == 6
    P = [-0.05, -0.1, -0.08, -0.1, -0.02, -0.04];
    Abar=[A, zeros(size(A,1),size(C,1));-C ,zeros(size(C,1),size(C,1))];
    Bbar=[B; zeros(size(C,1),2)];

    K = place(Abar, Bbar, P);
    K = vpa(K, 3);
    K_dyn = K(:,1:4);
    Kq = K(:, 5:6);
end

%% observer (full order)
FObsv_Poles = [-0.3, -0.4, -0.5, -0.6];
L = place(A', C', FObsv_Poles)';


%% LQR
Q = diag([500 500 0 0]);
R = diag([0.01 0.01]);
[K_lqr,S_lqr,e_lqr] = lqr(A,B,Q,R);
