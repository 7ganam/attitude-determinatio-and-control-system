function f=euler_wmap_fun(x,torq,d,jtrue)

% Pre-Allocate Space
f=zeros(7,1);torq=torq(:);d=d(:);

% Quaternion, Angular Velocity and Wheel Speed
q=[x(1);x(2);x(3);x(4)];
w=[x(5);x(6);x(7)];

% Integration Functions
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
   -w(2) w(1) 0];
om=[-wc w;-w' 0];
f(1:4,:)=0.5*om*q;
f(5:7,:)=inv(jtrue)*(-wc*jtrue*w+torq+d);
