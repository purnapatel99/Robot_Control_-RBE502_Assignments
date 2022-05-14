clear;
clc;

% Part b
syms th1 th2 dth1 dth2 ddth1 ddth2 u1 u2;
syms m1 m2 d1 d2 l1 l2 I1 I2 g T1 T2;

q = [th1; th2];
dq = [dth1; dth2];
ddq = [ddth1; ddth2];
u = [u1; u2];

% Manipulator form

M = [(m1*d1^2 + m2*d2^2 + 2*m2*cos(th2)*d2*l1 + m2*l1^2 + I1 + I2), (m2*d2^2 + l1*m2*cos(th2)*d2 + I2);
    (m2*d2^2 + l1*m2*cos(th2)*d2 + I2), (m2*d2^2 + I2)];

C = [-(2*d2*dth2*l1*m2*sin(th2)), -(d2*dth2*l1*m2*sin(th2));
    (d2*l1*m2*sin(th2)*dth1), 0];

G = [(- sin(th1)*(d1*g*m1 + g*l1*m2) - d2*g*m2*sin(th1 + th2));
    (- d2*g*m2*sin(th1 + th2))];

u = M*ddq+C*dq+G;

% Part c
% Feedback Linearization taking v = ddq

A = [0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
B = [0 0;
    0 0;
    1 0;
    0 1];

%lamda = [-1, -2, -2+1i, -2-1i];
% lamda = [-2, -5, -4+1i, -4-1i];
lamda = [-2, -5, -8+1i, -8-1i];

k = place(A,B,lamda)

syms des_th1 des_th2 des_dth1 des_dth2 vd1 vd2;

states = [q; dq];
des_states = [des_th1; des_th2; des_dth1; des_dth2];
vd = [vd1; vd2];

v = -k*(states-des_states)+vd;

% Final controller

tau = M*v+C*dq+G;

%Part e

%This function will run a visual simulation suimultaniously.
%In order to turn off the visual simulation, comment respective commands in the ode_RRbot.m file.

[t, y] = ode45(@ode_RRbot, [0,10], [deg2rad(200); deg2rad(125); 0; 0], [0; 0]);

b =  size(y);
T = [];

for a = 1: b(1)
    
    T(a,:)= -y(a,:)*(k');
    
end

des_th1 = (pi*t.^3)/500 - (3*pi*t.^2)/100 + pi;
des_th2 = (pi*t.^3)/1000 - (3*pi*t.^2)/200 + pi/2;
des_dth1 = (3*pi*t.^2)/500 - (3*pi*t)/50;
des_dth2 = (3*pi*t.^2)/1000 - (3*pi*t)/100;

% Plotting the output
figure;
plot(t,y(:,1));
hold on;
plot(t,des_th1);
xlabel('time (sec)');
ylabel('th1 (radian)');

figure;
plot(t,y(:,2));
hold on;
plot(t,des_th2);
xlabel('time (sec)');
ylabel('th2 (radian)');

figure;
plot(t,y(:,3));
hold on;
plot(t,des_dth1);
xlabel('time (sec)');
ylabel('dth1 (radian/sec)');

figure;
plot(t,y(:,4));
hold on;
plot(t,des_dth2);
xlabel('time (sec)');
ylabel('dth2 (radian/sec)');

figure;
plot(t,T(:,1));
xlabel('time (sec)');
ylabel('u1 (N.m)');

figure
plot(t,T(:,2));
xlabel('time (sec)');
ylabel('u2 (N.m)');





