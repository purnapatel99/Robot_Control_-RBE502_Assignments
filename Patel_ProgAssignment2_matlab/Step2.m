%Part b) c) d) e) and f)

clc
clear
m1=1; m2=1; l1=1; l2=1; d1=0.45; d2=0.45;
I1=0.084; I2=0.084; g=9.81;
syms th1 th2 dth1 dth2 ddth1 ddth2 T1 T2;


x = [th1; th2; dth1; dth2];
u = [T1; T2];

dX1=dth1;
dX2=dth2;

dX3=(I2*T1 - I2*T2 + T1*d2^2*m2 - T2*d2^2*m2 + d2^3*dth1^2*l1*m2^2*sin(th2)...
    + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1)...
    - T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2)...
    + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2)...
    + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1)...
    + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2...
    + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);

dX4=-(I2*T1 - I1*T2 - I2*T2 - T2*d1^2*m1 + T1*d2^2*m2 - T2*d2^2*m2 - T2*l1^2*m2...
    + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2)...
    - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1)...
    + I2*d1*g*m1*sin(th1) + T1*d2*l1*m2*cos(th2) - 2*T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1)...
    + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2)...
    + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2)...
    + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2)...
    + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1)...
    + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2)...
    + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(I1*I2 ...
    + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2 + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);

dX=[dX1;dX2;dX3;dX4];

% Part b

%leanirazation
A = jacobian (dX, x.');
B = jacobian (dX, u.');

% Part c

%subsituting equilibriun point (0,0)
A1 = subs(A,[th1, th2, dth1, dth2], [0,0,0,0]);
B1 = subs(B,[th1, th2, dth1, dth2], [0,0,0,0]);
A1=double(A1);
B1=double(B1);
%Stability
eigenA1=eig(A1)


%subsituting equilibriun point (pi,0)
A2 = subs(A,[th1, th2, dth1, dth2], [pi,0,0,0]);
B2 = subs(B,[th1, th2, dth1, dth2], [pi,0,0,0]);
A2=double(A2);
B2=double(B2);
%Stability
eigenA2=eig(A2)


%subsituting equilibriun point(0,pi)
A3 = subs(A,[th1, th2, dth1, dth2], [0,pi,0,0]);
B3 = subs(B,[th1, th2, dth1, dth2], [0,pi,0,0]);
A3=double(A3);
B3=double(B3);
%Stability
eigenA3=eig(A3)


%subsituting equilibriun point (pi,pi)
A4 = subs(A,[th1, th2, dth1, dth2], [pi,pi,0,0]);
B4 = subs(B,[th1, th2, dth1, dth2], [pi,pi,0,0]);
A4=double(A4);
B4=double(B4);
%Stability
eigenA4=eig(A4)

%Part d

%Controllability at upward configuration
C = ctrb(A1, B1)

rank = rank(ctrb(A1,B1))

%Part e

%lamda = [-5, -2, -2+1i, -2-1i];
lamda = [-1, -2, -2+1i, -2-1i];

k = place(A1,B1,lamda)

%Part f

%This function will run a visual simulation suimultaniously.
%In order to turn off the visual simulation, comment respective commands in the ode_RRbot.m file.

[t, y] = ode45(@ode_RRbot, [0,10], [pi/6; pi/4; 0; 0], [0; 0]);

b =  size(y);
T = [];

for a = 1: b(1)
    
    T(a,:)= -y(a,:)*(k');
    
end

% Plotting the output
figure;
plot(t,y(:,1));
xlabel('time (sec)');
ylabel('th1 (radian)');

figure;
plot(t,y(:,2));
xlabel('time (sec)');
ylabel('th2 (radian)');

figure;
plot(t,y(:,3));
xlabel('time (sec)');
ylabel('dth1 (radian/sec)');

figure;
plot(t,y(:,4));
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







