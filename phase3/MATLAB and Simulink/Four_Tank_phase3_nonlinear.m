%% Four Tank System_Phase 3 (Non-linear)
clc;
clear;
close all;

A = [730 730 730 730]; %% in cm^2
a = [2.10 2.14 2.2 2.3]; %% in cm^2
% Kc = 0.5; %% in V/cm
g = 981; %% in cm/s^2
y = [0.3 0.35];
k = [7.45 7.30];
 
syms h1 h2 h3 h4 h1_d h2_d h3_d h4_d v1 v2
h1_d = -a(1)/A(1)*sqrt(2*g*h1)+a(3)/A(1)*sqrt(2*g*h3)+y(1)*k(1)/A(1)*v1;
h2_d = -a(2)/A(2)*sqrt(2*g*h2)+a(4)/A(2)*sqrt(2*g*h4)+y(2)*k(2)/A(2)*v2;
h3_d = -a(3)/A(3)*sqrt(2*g*h3)+(1 - y(2))*k(2)/A(3)*v2;
h4_d = -a(4)/A(4)*sqrt(2*g*h4)+(1 - y(1))*k(1)/A(4)*v1;

F = [h1_d;h2_d;h3_d;h4_d];
h1_eq = 11.4;
h2_eq = 11.6;
equ1 = subs(h1_d,[h1,h2],[h1_eq h2_eq]);
equ2 = subs(h2_d,[h1,h2],[h1_eq h2_eq]);
equ3 = subs(h3_d,[h1,h2],[h1_eq h2_eq]);
equ4 = subs(h4_d,[h1,h2],[h1_eq h2_eq]);

answer = vpasolve([equ1==0;equ2==0;equ3==0; equ4==0],[v1,v2,h3,h4]);
v1_eq = real(double(answer.v1));
v2_eq = real(double(answer.v2));
% h1_eq = double(answer.h1)
% h2_eq = double(answer.h2)
h3_eq = real(double(answer.h3));
h4_eq = real(double(answer.h4));


Amat = jacobian(F,[h1,h2,h3,h4]);
Amat = subs(Amat,[h1,h2,h3,h4,v1,v2],[h1_eq,h2_eq,h3_eq,h4_eq,v1_eq,v2_eq]);
Amat = double(Amat);
Bmat = jacobian(F,[v1,v2]);
Bmat = subs(Bmat,[h1,h2,h3,h4,v1,v2],[h1_eq,h2_eq,h3_eq,h4_eq,v1_eq,v2_eq]);
Bmat = double(Bmat);
Cmat = [1 0 0 0;
   0 1 0 0];
Dmat = [0 0 ;0 0];