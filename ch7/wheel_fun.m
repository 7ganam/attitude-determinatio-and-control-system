function f=wheel_fun(x,jtrue,torq)

% Pre-Allocate Space
f=zeros(10,1);torq=torq(:);

% Quaternion, Angular Velocity and Wheel Momenta
q=[x(1);x(2);x(3);x(4)];
w=[x(5);x(6);x(7)];
wh=[x(8);x(9);x(10)];

% Integration Functions
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
   -w(2) w(1) 0];
om=[-wc w;-w' 0];
f(1:4,:)=0.5*om*q;
% f(5:7,:)=inv(jtrue-jw)*(-wc*(jtrue*w+jw*wh)+torq);
% f(8:10,:)=-inv(jw)*torq-f(5:7,:);
f(5:7,:)=inv(jtrue)*(-wc*jtrue*w+torq);
f(8:10,:)=-wc*wh-torq;