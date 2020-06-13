function [f,torq,hh,qm]=euler_paper_fun(x,jtrue,kp,kd,h,delta)

% Pre-Allocate Space
f=zeros(7,1);

% Quaternion, Angular Velocity and Wheel Speed
q=[x(1);x(2);x(3);x(4)];
w=[x(5);x(6);x(7)];

% Noisy Quaternion
m=0.2*rand(1);
ebar=randn(4,1);e=ebar/norm(ebar);
qm=(q+m*e)/norm(q+m*e);

if h*qm(4) >= -delta
 h=h;
else
 h=sign(qm(4));
end
hh=h;

% Control Torque 
torq=-kp*hh*qm(1:3)-kd*w;

% Integration Functions
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
   -w(2) w(1) 0];
om=[-wc w;-w' 0];
f(1:4,:)=0.5*om*q;
f(5:7,:)=inv(jtrue)*(-wc*jtrue*w+torq);
