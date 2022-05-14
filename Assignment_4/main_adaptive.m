 clear;
clc;

%Adaptive controller

%initial values
m1=0.75; m2=0.75; l1=1; l2=1; r1=0.45; r2=0.45;
I1=0.063; I2=0.063; g=9.81;

alpha_hat = [m2*l1^2 + m1*r1^2 + m2*r2^2 + I1 + I2;
        m2*l1*r2
        m2*r2^2 + I2
        m1*r1 + m2*l1
        m2*r2];


A = [0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
B = [0 0;
    0 0;
    1 0;
    0 1];


lamda = [-3, -3, -4, -4];

k = place(A,B,lamda);

Acl = A-B*k;
Q = eye(4);
P = lyap(Acl', Q);
gamma = eye(5).*0.1;


x0 = [deg2rad(200); deg2rad(125); 0; 0;0;0;alpha_hat(1);alpha_hat(2);alpha_hat(3);alpha_hat(4);alpha_hat(5)];

[t, y] = ode45(@(t,z) ode_adaptive_RRbot(t,z,k,P,gamma,B), [0,10], x0, [0; 0]);

b =  size(y);
T = [];

for a = 2: b(1)
    
    T(a,1)= (y(a,5)-y(a-1,5))/(t(a)-t(a-1));
    T(a,2)= (y(a,6)-y(a-1,6))/(t(a)-t(a-1));
    
    
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

figure
plot(t,y(:,7:11));
xlabel('time (sec)');
ylabel('alphahat');





