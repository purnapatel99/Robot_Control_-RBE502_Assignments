clc
clear

% Part a 

syms theta1 theta2 theta_dot1 theta_dot2 theta_ddot1 theta_ddot2 u1 u2;
syms m1 m2 r1 r2 l1 l2 I1 I2 g T1 T2;
q = [theta1; theta2];
dq = [theta_dot1; theta_dot2];
ddq = [theta_ddot1; theta_ddot2];
u = [u1; u2];
K = 0.5*(m1*r1^2 + I1 + m2*l1^2 + I2)*(theta_dot1^2) + 0.5*(m2*r2^2 + I2)*(theta_dot2^2) + (m2*r2^2 +I2)*theta_dot1*theta_dot2+ (m2*l1*r2)*((theta_dot1)^2 + theta_dot1*theta_dot2)*cos(theta2);
P = (m1*g*r1 + m2*g*l1)*cos(theta1) + (m2*g*r2)*cos(theta1 + theta2);
L = K - P;
DL_dq = jacobian(L, dq);
DL_q = jacobian(L, q);
dDL_dq = jacobian(DL_dq, [q; dq])*[dq; ddq];

u = [dDL_dq(1) - DL_q(1); dDL_dq(2) - DL_q(2)];

%Part b

X = sym ('X', [4,1]);
X(1) = theta1;
X(2) = theta2;
X(3) = theta_dot1;
X(4) = theta_dot2;

eq1 = u(1) - T1;
eq2 = u(2) - T2;

sol =  solve([eq1==0, eq2==0], [theta_ddot1, theta_ddot2]);
sol.theta_ddot1;
sol.theta_ddot2;

%Part c

[t, y] = ode45(@ode_RRbot, [0,10], [pi/6; pi/4; 0; 0]);
figure;
plot(t,y);
xlabel('time (sec)');
ylabel('angle (radian)');