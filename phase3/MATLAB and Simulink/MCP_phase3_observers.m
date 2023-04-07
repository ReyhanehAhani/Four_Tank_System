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

F = [h1d; h2d; h3d; h4d];

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


%% observer (full order)
FObsv_Poles = [-0.3, -0.4, -0.5, -0.6];
L = place(A', C', FObsv_Poles)';

%% observer (reduced order)

RObsv_Poles = [-0.4, -.45];   % We can't choose -0.5 cause it's exactly the eigvalue for A matrix
F_mat = diag(RObsv_Poles);
L_mat = [1, 0;0, -1];
T_mat = (lyap(A', -F_mat', -C'*L_mat'))';

