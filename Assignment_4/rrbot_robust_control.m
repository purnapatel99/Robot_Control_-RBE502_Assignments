%Part f
clear; close; clc;

% ROS Setup
rosinit;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');
tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;
send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
% Initilizating Paremeters
t = 0;

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
p = [10 0; 0 10]; %rho
phi = 0.05;

m1=0.75; m2=0.75; l1=1; l2=1; d1=0.45; d2=0.45;
I1=0.063; I2=0.063; g=9.81;

x = [];
u = [];
t1 = [];
n = 1;

while(t < 10)
    t = toc;
    
    % Trajectory
    
    des_state = [(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
            (pi*t^3)/1000 - (3*pi*t^2)/200 + pi/2;
            (3*pi*t^2)/500 - (3*pi*t)/50;
            (3*pi*t^2)/1000 - (3*pi*t)/100];
        
    vd = [(3*pi*t)/250 - (3*pi)/50;
        (3*pi*t)/500 - (3*pi)/100];

    % read the joint states
    jointData = receive(JointStates);
    
    % inspect the "jointData" variableinMATLAB to get familiar with itsstructure
    
    state = [jointData.Position;jointData.Velocity];
    
    th1 = state(1);
    th2 = state(2);
    dth1 = state(3);
    dth2 = state(4);
    
    e = state - des_state;
    
    % Robust Controller
    
    if norm(e'*P*B)> phi
        vr = -((e'*P*B)/norm(e'*P*B))*p;
    else
        vr = -(e'*P*B)*p/phi;
    end

    v = vd -k*(e)+vr';

    M = [(m1*d1^2 + m2*d2^2 + 2*m2*cos(th2)*d2*l1 + m2*l1^2 + I1 + I2), (m2*d2^2 + l1*m2*cos(th2)*d2 + I2);
        (m2*d2^2 + l1*m2*cos(th2)*d2 + I2), (m2*d2^2 + I2)];

    C = [-(2*d2*dth2*l1*m2*sin(th2)), -(d2*dth2*l1*m2*sin(th2));
        (d2*l1*m2*sin(th2)*dth1), 0];

    G = [(- sin(th1)*(d1*g*m1 + g*l1*m2) - d2*g*m2*sin(th1 + th2));
        (- d2*g*m2*sin(th1 + th2))];

    T = M*v+C*[dth1; dth2]+G;
    
    tau1.Data = T(1);
    tau2.Data = T(2);
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    % Sampelling data
    
    x(:,n) = [jointData.Position;jointData.Velocity];
    u(:,n) = [jointData.Effort];
    t1(:,n) = t;
    
    n = n+1;
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

des_th1 = (pi*t1.^3)/500 - (3*pi*t1.^2)/100 + pi;
des_th2 = (pi*t1.^3)/1000 - (3*pi*t1.^2)/200 + pi/2;
des_dth1 = (3*pi*t1.^2)/500 - (3*pi*t1)/50;
des_dth2 = (3*pi*t1.^2)/1000 - (3*pi*t1)/100;

%Plotting the output
figure;
plot(t1,x(1,:));
hold on;
plot(t1,des_th1);
xlabel('time (sec)');
ylabel('th1 (radian)');

figure;
plot(t1,x(2,:));
hold on;
plot(t1,des_th2);
xlabel('time (sec)');
ylabel('th2 (radian)');

figure;
plot(t1,x(3,:));
hold on;
plot(t1,des_dth1);
xlabel('time (sec)');
ylabel('dth1 (radian/sec)');

figure;
plot(t1,x(4,:));
hold on;
plot(t1,des_dth2);
xlabel('time (sec)');
ylabel('dth2 (radian/sec)');

figure;
plot(t1,u(1,:));
xlabel('time (sec)');
ylabel('u1 (N.m)');

figure
plot(t1,u(2,:));
xlabel('time (sec)');
ylabel('u2 (N.m)');

% disconnect from roscore
rosshutdown;







