% Stability of four tank system
clc;
clear;
close all;

Q = eye(4);    % positive definite matrix
A =[ -0.01887,        0,  0.02899,        0;
         0, -0.01906,        0,  0.03489;
        0,        0, -0.02899,        0;
        0,        0,        0, -0.03489];

P = vpa(lyap(A, Q), 4);
D1 = det(P(1,1));
D2 = det(P(1:2,1:2));
D3 = det(P(1:3,1:3));
D4 = det(P);
if D1>0 && D2>0 && D3>0 && D4>0
    disp('System is stable')
elseif D1<0 && D2>0 && D3<0 && D4>0
    disp('System is not stable')
end
