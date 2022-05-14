function dz = ode_RRbot(t, z)

m1=1; m2=1; l1=1; l2=1; r1=0.45; r2=0.45;
I1=0.084; I2=0.084; g=9.81;

dz = zeros (4,1);
z = num2cell(z);
[theta1, theta2, theta_dot1, theta_dot2] = deal(z{:});

if abs(theta1)>2*pi
    theta1 = mod(theta1, 2*pi);
end
if abs(theta2)>2*pi
    theta2 = mod(theta2, 2*pi);
end

T1=0;
T2=0;   
dz(1)=theta_dot1;
dz(2)=theta_dot2;
dz(3)=(I2*T1 - I2*T2 + T1*m2*r2^2 - T2*m2*r2^2 + l1*m2^2*r2^3*theta_dot1^2*sin(theta2) + l1*m2^2*r2^3*theta_dot2^2*sin(theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) - T2*l1*m2*r2*cos(theta2) + 2*l1*m2^2*r2^3*theta_dot1*theta_dot2*sin(theta2) + l1^2*m2^2*r2^2*theta_dot1^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + I2*l1*m2*r2*theta_dot1^2*sin(theta2) + I2*l1*m2*r2*theta_dot2^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*I2*l1*m2*r2*theta_dot1*theta_dot2*sin(theta2))/(I1*I2 - m2^2*r2^4 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 - I2*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2 + m1*m2*r1^2*r2^2);
dz(4)=-(I2*T1 - I1*T2 - I2*T2 - T2*l1^2*m2 - T2*m1*r1^2 + T1*m2*r2^2 + g*m2^2*r2^3*sin(theta1 + theta2) + l1^3*m2^2*r2*theta_dot1^2*sin(theta2) + l1*m2^2*r2^3*theta_dot2^2*sin(theta2) - g*l1^2*m2^2*r2*sin(theta1 + theta2) - I1*g*m2*r2*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1) + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + T1*l1*m2*r2*cos(theta2) - 2*T2*l1*m2*r2*cos(theta2) + 2*l1*m2^2*r2^3*theta_dot1*theta_dot2*sin(theta2) + 2*l1^2*m2^2*r2^2*theta_dot1^2*cos(theta2)*sin(theta2) + l1^2*m2^2*r2^2*theta_dot2^2*cos(theta2)*sin(theta2) - g*l1*m2^2*r2^2*sin(theta1 + theta2)*cos(theta2) + g*l1^2*m2^2*r2*cos(theta2)*sin(theta1) - g*m1*m2*r1^2*r2*sin(theta1 + theta2) + I1*l1*m2*r2*theta_dot1^2*sin(theta2) + I2*l1*m2*r2*theta_dot1^2*sin(theta2) + I2*l1*m2*r2*theta_dot2^2*sin(theta2) + g*m1*m2*r1*r2^2*sin(theta1) + 2*l1^2*m2^2*r2^2*theta_dot1*theta_dot2*cos(theta2)*sin(theta2) + l1*m1*m2*r1^2*r2*theta_dot1^2*sin(theta2) + 2*I2*l1*m2*r2*theta_dot1*theta_dot2*sin(theta2) + g*l1*m1*m2*r1*r2*cos(theta2)*sin(theta1))/(I1*I2 - m2^2*r2^4 + l1^2*m2^2*r2^2 + I2*l1^2*m2 + I2*m1*r1^2 + I1*m2*r2^2 - I2*m2*r2^2 - l1^2*m2^2*r2^2*cos(theta2)^2 + m1*m2*r1^2*r2^2);

end
