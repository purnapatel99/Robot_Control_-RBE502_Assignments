%Part e
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
req.JointPositions = [deg2rad(30), deg2rad(45)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;

% k for lamda [-5, -2, -2+1i, -2-1i]
%k = [41.0473   13.2298   15.9180    6.0257; 11.3837    7.5267    4.9180    2.1538];

% k for lamda [-1, -2, -2+1i, -2-1i]
k = [24.9780    5.0173    7.8833    1.9194; 6.1250    4.8392    2.2887    0.8101];

x = [];
u = [];
t1 = [];
n = 1;

while(t < 10)
    t = toc;
    
    % read the joint states
    jointData = receive(JointStates);
    
    % inspect the "jointData" variableinMATLAB to get familiar with itsstructure
    
    % design your state feedback controllerinthe following
    tau1.Data = -k(1,:)*[jointData.Position;jointData.Velocity];
    tau2.Data = -k(2,:)*[jointData.Position;jointData.Velocity];
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);
    
    % you can sample data here to plot at the end
    
    x(:,n) = [jointData.Position;jointData.Velocity];
    u(:,n) = [jointData.Effort];
    t1(:,n) = t;
    
    n = n+1;
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

%Plotting the output
figure;
plot(t1,x(1,:));
xlabel('time (sec)');
ylabel('th1 (radian)');

figure;
plot(t1,x(2,:));
xlabel('time (sec)');
ylabel('th2 (radian)');

figure;
plot(t1,x(3,:));
xlabel('time (sec)');
ylabel('dth1 (radian/sec)');

figure;
plot(t1,x(4,:));
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







