clear;
clc;

% Part a

% Defining initial and final Conditions
t0 = 0; tf = 10;
th1_0 = deg2rad(180);th1_f = 0;
th2_0 = deg2rad(90); th2_f = 0;
dth1_0 = 0; dth1_f = 0;
dth2_0 = 0; dth2_f = 0;
syms t;

A = [1 t0 t0^2 t0^3;
    0 1 2*t0 3*t0^2;
    1 tf tf^2 tf^3;
    0 1 2*tf 3*tf^2];

B = [th1_0, th2_0; dth1_0, dth2_0; th1_f, th2_f; dth1_f, dth2_f];

a = inv(A)*B;

th1 = a(1,1) + a(2,1)*t + a(3,1)*t^2 + a(4,1)*t^3
th2 = a(1,2) + a(2,2)*t + a(3,2)*t^2 + a(4,2)*t^3
dth1 = jacobian(th1, t)
dth2 = jacobian(th2, t)
ddth1 = jacobian(dth1, t)
ddth2 = jacobian(dth2, t)