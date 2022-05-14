clear;
clc;

%Robust controller

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
Q = eye(4).*5;
P = lyap(Acl', Q);
p = [10 0; 0 10]; %rho
phi = 0.025;


[t, y] = ode45(@(t,z) ode_robust_RRbot(t,z,k,P,p,B,phi), [0,10], [deg2rad(200); deg2rad(125); 0; 0;0;0], [0; 0]);


%Calculating u1 and u2
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





