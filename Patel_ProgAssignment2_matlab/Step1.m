%Part a

clc
clear
m1=1; m2=1; l1=1; l2=1; d1=0.45; d2=0.45;
I1=0.084; I2=0.084; g=9.81;
syms th1 th2 dth1 dth2 ddth1 ddth2 T1 T2;
sol = zeros(4, 2);


x = [th1, th2, dth1, dth2];
u = [T1, T2];
syms dX;


dX(1)=dth1;

dX(2)=dth2;

dX(3)=(I2*T1 - I2*T2 + T1*d2^2*m2 - T2*d2^2*m2 + d2^3*dth1^2*l1*m2^2*sin(th2)...
    + d2^3*dth2^2*l1*m2^2*sin(th2) + d2^2*g*l1*m2^2*sin(th1) + I2*d1*g*m1*sin(th1)...
    - T2*d2*l1*m2*cos(th2) + I2*g*l1*m2*sin(th1) + 2*d2^3*dth1*dth2*l1*m2^2*sin(th2)...
    + d2^2*dth1^2*l1^2*m2^2*cos(th2)*sin(th2) - d2^2*g*l1*m2^2*sin(th1 + th2)*cos(th2)...
    + I2*d2*dth1^2*l1*m2*sin(th2) + I2*d2*dth2^2*l1*m2*sin(th2) + d1*d2^2*g*m1*m2*sin(th1)...
    + 2*I2*d2*dth1*dth2*l1*m2*sin(th2))/(I1*I2 + d2^2*l1^2*m2^2 + I2*d1^2*m1 + I1*d2^2*m2...
    + I2*l1^2*m2 + d1^2*d2^2*m1*m2 - d2^2*l1^2*m2^2*cos(th2)^2);

dX(4)=-(I2*T1 - I1*T2 - I2*T2 - T2*d1^2*m1 + T1*d2^2*m2 - T2*d2^2*m2 - T2*l1^2*m2...
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

dX(3)=subs(dX(3),[dth1,dth2,T1,T2],[0,0,0,0]);
dX(4)=subs(dX(4),[dth1,dth2,T1,T2],[0,0,0,0]);

dX(1)=x(1)>=0;
dX(2)=x(2)>=0;

% The solve function gives only the first equilibriun available
%please use newer version of matlab, because when in the older version
%solve function cannot take 4 conditions

[sol(1,1), sol(1,2)]=solve(dX(3)==0, dX(4)==0, dX(1), dX(2));

% For loop is implements different conditions in order to find all the
% equilibrium points

for a = 1:2
    
    for b = 1:2
        dX(1)=x(1)>=sol(1,1);
        dX(2)=x(2)>=sol(1,2);
        dX(a)=x(a)>sol(1,a);
        dX(b)=x(b)>sol(1,b);
        [sol(a+b,1), sol(a+b,2)]=solve(dX(3)==0, dX(4)==0, dX(1), dX(2));
        
    end
    
end
    
sol    
    
    



