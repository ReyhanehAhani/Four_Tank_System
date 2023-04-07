%% Four Tank System_Phase 1
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
taw =[2 2.1];             % pump response time constants in s

%% Non-linear differential equations
syms h1 h2 h3 h4 h1_d h2_d h3_d h4_d v1 v2 u1 u2
h1d = -a(1)/Area(1)*sqrt(2*g*h1)+a(3)/Area(1)*sqrt(2*g*h3)+y(1)*k(1)/Area(1)*v1+0.00137;
h2d = -a(2)/Area(2)*sqrt(2*g*h2)+a(4)/Area(2)*sqrt(2*g*h4)+y(2)*k(2)/Area(2)*v2+0.00137;
h3d = -a(3)/Area(3)*sqrt(2*g*h3)+(1 - y(2))*(k(2)/Area(3))*v2;
h4d = -a(4)/Area(4)*sqrt(2*g*h4)+(1 - y(1))*(k(1)/Area(4))*v1;
v1d = -v1/taw(1)+u1/taw(1);
v2d = -v2/taw(2)+u2/taw(2);

%% Linearization and Steady-State model
F = [h1d;h3d;v2d;v1d;h2d;h4d];

Amat = jacobian(F,[h1,h3,v2,v1,h2,h4]);
Amat = subs(Amat,[h1,h3,v2,v1,h2,h4,u1,u2],[h1_eq,h3_eq,v2_eq,v1_eq,h2_eq,h4_eq,0,0]);
A = double(Amat);

Bmat = jacobian(F,[u1,u2]);
Bmat = subs(Bmat,[h1,h2,h3,h4,v1,v2,u1,u2],[h1_eq,h2_eq,h3_eq,h4_eq,v1_eq,v2_eq,0,0]);
B = double(Bmat);

C = [1 0 0 0 0 0;
     0 0 0 0 1 0];
 
D = [0 0 ;0 0];
sys = ss(A,B,C,D)

%% Transfer functions 
G = tf(sys)

%% Pole/Zero map of transfer functions
figure('Name','Pole/Zero map of system','NumberTitle','off')
pzmap(G)
Z = tzero(sys)        % Invariant zeros of linear system

figure('Name','Pole/Zero map','NumberTitle','off')
k = 1;
for i = 1:size(G,1)
    for j = 1:size(G,2)
        subplot(2,2,k)
        pzmap(G(i,j))
        title(['G(',num2str(i),',',num2str(j),')'])
        k = k+1;
    end
    
end
%% Root locus of tranfer functions
figure('Name','Root Locus','NumberTitle','off')
subplot(2,2,1)
rlocus(G(1,1))
title("G(1,1)")

subplot(2,2,2)
rlocus(G(1,2))
title("G(1,2)")

subplot(2,2,3)
rlocus(G(2,1))
title("G(2,1)")

subplot(2,2,4)
rlocus(G(2,2))
title("G(2,2)")

%% State Transition Matrix
syms s
phi = inv((s*eye(6)-A));
phi = simplify(phi);
phi_s = vpa(phi, 3)

syms t
phi_t = simplify(expm(A*t));
phi_t = vpa(phi_t,3)

%% Step Response of system
figure('Name','Step Response of system','NumberTitle','off')
step(sys)
grid on

%% Initial condition response of system
x0 = rand(6,1)
figure('Name','Initial condition response of system','NumberTitle','off')
initial(sys,x0)
grid on

%% controllability & observability of the system
Co = ctrb(A, B)
Rc = rank(Co)
Ob = obsv(A, C)
Ro = rank(Ob)
% Since both ranks are 6, this system is controllable & observable

%% Jordan form (J) & the Similarity transform (V)
[V,J] = jordan(A)
P = inv(V)

%% Initial condition due to non-excite specific frequency
syms x01 x02 x03 x04 x05 x06
x00 = [x01 x02 x03 x04 x05 x06].';
cc = simplify(V \ x00);             % V\x00 = inv(V)*x00 = P*x00
cc = vpa(cc,3)

%% Kalman decomposition
[kalman_sys, U] = minreal(sys);
A_minreal = kalman_sys.A
B_minreal = kalman_sys.B
C_minreal = kalman_sys.C


