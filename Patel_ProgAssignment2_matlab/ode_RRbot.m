function dz = ode_RRbot(t, z)

m1=1; m2=1; l1=1; l2=1; d1=0.45; d2=0.45;
I1=0.084; I2=0.084; g=9.81;

dz = zeros (4,1);
z = num2cell(z);

[th1, th2, dth1, dth2] = deal(z{:});

if abs(th1)>2*pi
    th1 = mod(th1, 2*pi);
end
if abs(th2)>2*pi
    th2 = mod(th2, 2*pi);
end

% k for lamda [-1, -2, -2+1i, -2-1i]
k = [24.9780    5.0173    7.8833    1.9194;
    6.1250    4.8392    2.2887    0.8101];


% k for lamda [-5, -2, -2+1i, -2-1i]
%k = [41.0473   13.2298   15.9180    6.0257;
%    11.3837    7.5267    4.9180    2.1538];

T = -k*[th1; th2; dth1; dth2];

T1=T(1,1);
T2=T(2,1);
dz(1)=dth1;
dz(2)=dth2;
dz(3)=(I2*T1 - I2*T2 + T1*d2^2*m2 - T2*d2^2*m2 + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) - T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2 + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);
dz(4)=-(I2*T1 - I1*T2 - I2*T2 - T2*d1^2*m1 + T1*d2^2*m2 - T2*d2^2*m2 - T2*l1^2*m2 + d2*dth1^2*l1^3*m2^2*sin(th2) + d2^3*dth1^2*l1*m2^2*sin(th2) + d2^3*dth2^2*l1*m2^2*sin(th2) - d2*g*l1^2*m2^2*sin(th1 + th2) - I1*d2*g*m2*sin(th1 + th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1) + T1*d2*l1*m2*cos(th2) - 2*T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2) + 2*d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) + d2^2*dth2^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2) + d2*g*l1^2*m2^2*cos(th2)*sin(th1) - d1^2*d2*g*m1*m2*sin(th1 + th2) + I1*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1) + 2*d2^2*dth1*dth2*l1^2*m2^2*cos(th2)*sin(th2) + d1^2*d2*dth1^2*l1*m1*m2*sin(th2) + 2*I2*d2*dth1*dth2*l1*m2*sin(th2) + d1*d2*g*l1*m1*m2*cos(th2)*sin(th1))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2 + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);


%For running visual simulation. Comment below this to stop visual simulation.
x1 = l1*sin(th1);
y1 = l1*cos(th1);
x2 = l1*sin(th1) + l2*sin(th1+th2);
y2 = l1*cos(th1) + l2*cos(th1+th2);
plot(2,2);
plot(-2,-2);
hold on;
plot(2,2);
plot(-2,-2);
plot(0,0,'r.');
plot([0 x1], [0 y1]);
plot([x1 x2], [y1 y2]);
plot(x1, y1, 'r.');
plot(x2, y2, 'r.');
hold off
pause(1/1000);

end
