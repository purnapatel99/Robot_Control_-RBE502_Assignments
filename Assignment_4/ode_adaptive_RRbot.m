function [dz] = ode_adaptive_RRbot(t, z, k, P, gamma, B)
dz = zeros (11,1);
z = num2cell(z);

[th1, th2, dth1, dth2, T1, T2, alpha_hat(1,1), alpha_hat(2,1), alpha_hat(3,1), alpha_hat(4,1), alpha_hat(5,1)] = deal(z{:});

if abs(th1)>2*pi
    th1 = mod(th1, 2*pi);
end
if abs(th2)>2*pi
    th2 = mod(th2, 2*pi);
end

%Part d

state = [th1; th2; dth1; dth2];

%Desired States from trajectory

des_state = [(pi*t^3)/500 - (3*pi*t^2)/100 + pi;
            (pi*t.^3)/1000 - (3*pi*t.^2)/200 + pi/2;
            (3*pi*t^2)/500 - (3*pi*t)/50;
            (3*pi*t^2)/1000 - (3*pi*t)/100];
        
vd = [(3*pi*t)/250 - (3*pi)/50;
    (3*pi*t)/500 - (3*pi)/100];

e = state - des_state;


% k for lamda = [-3,-3,-4,-4]

v = vd -k*(e);
g = 9.81;
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

%%

%Defining the original values for the simulation
m1=1; m2=1; l1=1; l2=1; d1=0.45; d2=0.45;
I1=0.084; I2=0.084; g=9.81;

T1=T(1);
T2=T(2);
dz(1)=dth1;
dz(2)=dth2;
ddth1=(I2*T1 - I2*T2 + T1*d2^2*m2 - T2*d2^2*m2 + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) - T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2 + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);
ddth2=-(I2*T1 - I1*T2 - I2*T2 - T2*d1^2*m1 + T1*d2^2*m2 - T2*d2^2*m2 - T2*l1^2*m2 + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + T1*d2*l1*m2*cos(th2) - 2*T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2 + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);
dz(3)= ddth1;
dz(4)= ddth2;
dz(5)= T1;
dz(6)= T2;

%%
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

dz(7) = d_alpha_hat(1);
dz(8) = d_alpha_hat(2);
dz(9) = d_alpha_hat(3);
dz(10) = d_alpha_hat(4);
dz(11) = d_alpha_hat(5);
%ans = t


end
