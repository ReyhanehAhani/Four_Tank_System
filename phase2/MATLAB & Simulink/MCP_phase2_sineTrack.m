clc;
clear;
close all;
%% Sine tracking 
% input 1 to output 2
num = 0.0001246;
den = [1 0.554 0.02764 0.0003325];
[A, B, C, D] = tf2ss(num, den);

factors = factor(poly2sym(den));
OL_Poles = roots(den);

%% without disturbance
num_1 = num;
den_1 = conv(den, [1 0 1]);
[A1, B1, C1, D1] = tf2ss(num_1, den_1);

K1 = place(A1,B1, [-0.025, -0.05, -0.25, -0.9, -1]); 
L1 = place(A1', C1', [-5, -3, -4, -6, -4.5])';

% K1 controller
[b1, a1] = ss2tf(A1-B1*K1-L1*C1,L1,K1,0);

%% with disturbance
num_2 = num;
% Considering r_d & w dynamics (1/(s(s^2+1))
den_2 = conv(den, [1 0 1 0]);  
[A2, B2, C2, D2] = tf2ss(num_2, den_2);
 
K2 = place(A2, B2, [-0.025, -0.05, -0.25, -0.9, -1, -1.1]); 
L2 = place(A2', C2', [-3,-4,-5,-6,-7,-8])';
[b2, a2] = ss2tf(A2-B2*K2-L2*C2, L2, K2, 0);