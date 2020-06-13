% This program simulates the second example of the following paper:
% Mayhew, C.G., Sanfelice, R.G., and Teel, A.R., "Quaternion-Based Hybrid
% Control for Robust Global Attitude Tracking," IEEE Transactions on 
% Automatic Control, Vol. 56, No. 11, Nov. 2011, pp. 2555-2566.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.7

% Written by John L. Crassidis 9/13

% Other Required Routines: euler_paper_fun.m

% v Value
v=[1;2;3]/norm([1;2;3]);

% Inertia 
in=10*diag(v);

% Time
dt=0.001;tf=70;
t=[0:dt:tf]';
m=length(t);

% Pre-Allocate Space
q=zeros(m,4);
w=zeros(m,3);
u=zeros(m,3);
xa=zeros(m,7);
i1000=0;

% Initial Quaternion and Rate
q(1,:)=[sqrt(1-0.2^2)*v' -0.2];
w(1,:)=0.5*v';
xa(1,:)=[q(1,:) w(1,:)];

% Gains
kp=1;kd=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discontinuous
h_store=zeros(m,1);
delta=0;
h=1;

% Get Torque and h at Initial Time
[ff,torq,hh]=euler_paper_fun(xa(1,:),in,kp,kd,h,delta);
u(1,:)=torq(:)';
h_store(1)=hh;
torq_norm_dis=zeros(m,1);

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

[f1,torq,hh]=euler_paper_fun(xa(i,:),in,kp,kd,h,delta);h=hh;
[f2,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f1'*dt,in,kp,kd,h,delta);h=hh;
[f3,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f2'*dt,in,kp,kd,h,delta);h=hh;
[f4,torq,hh]=euler_paper_fun(xa(i,:)+f3'*dt,in,kp,kd,h,delta);h=hh;
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Get Torque and h
[ff,torq,hh]=euler_paper_fun(xa(i+1,:),in,kp,kd,h,delta);
u(i+1,:)=torq(:)';
h_store(i+1)=hh;

% Actual Quaterion and Angular Velocity 
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Torque Norm Integral
torq_norm_dis(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));

end

h_store_q_dis=h_store.*q(:,4);
theta_dis=2*acos(abs(q(:,4)))*180/pi;
w_dis=w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unwinding (delta > 1) 
h_store=zeros(m,1);
delta=2;
h=1;

% Get Torque and h at Initial Time
[ff,torq,hh]=euler_paper_fun(xa(1,:),in,kp,kd,h,delta);
u(1,:)=torq(:)';
h_store(1)=hh;
torq_norm_unwind=zeros(m,1);

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

[f1,torq,hh]=euler_paper_fun(xa(i,:),in,kp,kd,h,delta);h=hh;
[f2,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f1'*dt,in,kp,kd,h,delta);h=hh;
[f3,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f2'*dt,in,kp,kd,h,delta);h=hh;
[f4,torq,hh]=euler_paper_fun(xa(i,:)+f3'*dt,in,kp,kd,h,delta);h=hh;
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Get Torque and h
[ff,torq,hh]=euler_paper_fun(xa(i+1,:),in,kp,kd,h,delta);
u(i+1,:)=torq(:)';
h_store(i+1)=hh;

% Actual Quaterion and Angular Velocity 
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Torque Norm Integral
torq_norm_unwind(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));

end

h_store_q_unwind=h_store.*q(:,4);
theta_unwind=2*acos(abs(q(:,4)))*180/pi;
w_unwind=w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid 
h_store=zeros(m,1);
delta=0.4;
h=1;
torq_norm=zeros(m,1);

% Get Torque and h at Initial Time
[ff,torq,hh]=euler_paper_fun(xa(1,:),in,kp,kd,h,delta);
u(1,:)=torq(:)';
h_store(1)=hh;

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

[f1,torq,hh]=euler_paper_fun(xa(i,:),in,kp,kd,h,delta);h=hh;
[f2,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f1'*dt,in,kp,kd,h,delta);h=hh;
[f3,torq,hh]=euler_paper_fun(xa(i,:)+0.5*f2'*dt,in,kp,kd,h,delta);h=hh;
[f4,torq,hh]=euler_paper_fun(xa(i,:)+f3'*dt,in,kp,kd,h,delta);h=hh;
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Get Torque and h
[ff,torq,hh]=euler_paper_fun(xa(i+1,:),in,kp,kd,h,delta);
u(i+1,:)=torq(:)';
h_store(i+1)=hh;

% Actual Quaterion and Angular Velocity 
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Torque Norm Integral
torq_norm(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));

end

h_store_q=h_store.*q(:,4);
theta=2*acos(abs(q(:,4)))*180/pi;

plot(t,h_store_q,t,h_store_q_dis,'--',t,h_store_q_unwind,'-.')
set(gca,'fontsize',12)
legend('Hybrid','Discontinuous','Unwinding','Location','SouthEast')
axis([0 15 -1 1])
ylabel('aaa')
xlabel('Time (Sec)')

pause

plot(t,theta,t,theta_dis,'--',t,theta_unwind,'-.')
set(gca,'fontsize',12)
legend('Hybrid','Discontinuous','Unwinding')
axis([0 70 0 180])
ylabel('Theta')
xlabel('Time (Sec)')

pause

plot(t,(w(:,1).^2+w(:,2).^2+w(:,3).^2).^(0.5),t,(w_dis(:,1).^2+w_dis(:,2).^2+w_dis(:,3).^2).^(0.5),'--',t,(w_unwind(:,1).^2+w_unwind(:,2).^2+w_unwind(:,3).^2).^(0.5),'-.')
set(gca,'fontsize',12)
legend('Hybrid','Discontinuous','Unwinding')
axis([0 70 0 0.6])
set(gca,'ytick',[0 0.1 0.2 0.3 0.4 0.5 0.6])
ylabel('Angular Velocity')
xlabel('Time (Sec)')

pause
plot(t,torq_norm,t,torq_norm_dis,'--',t,torq_norm_unwind,'-.')
set(gca,'fontsize',12)
legend('Hybrid','Discontinuous','Unwinding','Location','SouthEast')
axis([0 70 0 3.5])
ylabel('Torque Norm')
xlabel('Time (Sec)')