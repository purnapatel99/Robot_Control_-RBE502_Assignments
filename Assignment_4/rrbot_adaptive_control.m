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
t = 0;
t_prev = 0;
dth1_prev = 0;
dth2_prev = 0;

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
gamma = eye(5).*9;

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
    ddth1 = 1000*(dth1 - dth1_prev)/(1000*(t - t_prev));
    ddth2 = 1000*(dth2 - dth2_prev)/(1000*(t - t_prev));
    %ans = t - t_prev
    e = state - des_state;
    
    % Adaptive controller
    v = vd -k*(e);

    Y_out = [v(1), ...
        cos(th2)*(2*v(1) + v(2)) - 2*sin(th2)*dth1*dth2 - sin(th2)*dth2^2, ...
        v(2), ...
        -sin(th1)*g, ...
        -sin(th1 + th2)*g; ...
        0, ...
        sin(th2)*dth1^2 + cos(th2)*v(1), ...
        v(1) + v(2), ...
        0, ...
        -sin(th1+th2)*g];

    T = Y_out*alpha_hat;
    
    %Adaption law
    
    Y = [ddth1, ...
        cos(th2)*(2*ddth1 + ddth2) - 2*sin(th2)*dth1*dth2 - sin(th2)*dth2^2, ...
        ddth2, ...
        -sin(th1)*g, ...
        -sin(th1 + th2)*g; ...
        0, ...
        sin(th2)*dth1^2 + cos(th2)*ddth1, ...
        ddth1 + ddth2, ...
        0, ...
        -sin(th1+th2)*g];

    M1 = [1, 2*cos(th2), 0, 0, 0;
          0, cos(th2), 1, 0, 0];
    M2 = [0, cos(th2), 1, 0, 0;
          0, 0, 1, 0, 0];
    M_hat = [M1*alpha_hat, M2*alpha_hat];
    phi = M_hat\Y;
  
    d_alpha_hat = -gamma\(phi'*B'*P*e);
    
    %Integrating d_aplha_hat and updating alpha  
    alpha_hat = alpha_hat + d_alpha_hat*(t-t_prev);
    
    %Uncomment following lines to use ode45 for integrating
    %[tm,y] = ode45(@(tm,y) d_alpha_hat, [t_prev t], alpha_hat);
    %b = size(y);
    %alpha_hat = y(b(1), :)';
    
    tau1.Data = T(1);
    tau2.Data = T(2);
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    
    
    %Sampelling data
    
    x(:,n) = [jointData.Position;jointData.Velocity];
    u(:,n) = [jointData.Effort];
    t1(:,n) = t;
    t_prev = t;
    dth1_prev = dth1;
    dth2_prev = dth2;
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







