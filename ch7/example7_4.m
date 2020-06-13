% This example shows performs a reorientation manuever for the spacecraft
% regulation case involving thrusters.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.4

% Written by John L. Crassidis 9/13

% Other Required Routines: euler_quat_reg_thrust_fun.m, quat_err.m, schmidt.m

% Time and Inertia Matrix
dt=0.1;tf=300;
t=[0:dt:tf]';m=length(t);
j=diag([10000 9000 12000]);invj=inv(j);

% Initial and Desired Quaternion, and Initial Rate
qd=[0;0;0;1];
q0=[0.685;0.695;0.153;0.153];q0=q0/norm(q0);
w0=[0.53;0.53;0.053]*pi/180;
x=zeros(m,7);x(1,:)=[q0' w0'];

% Gains
kp=50;kd=500;

% Desired Quaternion Matrix
quat_mult_mat=[qd(4) qd(3) -qd(2) -qd(1)
              -qd(3) qd(4) qd(1) -qd(2)
               qd(2) -qd(1) qd(4) -qd(3)
               qd(1) qd(2) qd(3) qd(4)];

torque=zeros(m,3);

% PWPF Parameters
kpm=275;
km=5;
tau=0.5;
uon=10;
uoff=0.6*uon;
umax=20;

% Set Initial Torque to Zero
u=zeros(3,1);ufilt=u;

% Convert Filter to Discrete Time
num=km;den=[tau 1];
[af,bf,cf,df]=tf2ss(num,den);
[phi,gam]=c2d(af,bf,dt);

 % Main Loop
for i = 1:m-1
% Error Quaternion and Feedback Error Signal
 qerr=quat_mult_mat*x(i,1:4)';
 u_err=(-kp*sign(qerr(4))*qerr(1:3)-kd*(1+qerr(1:3)'*qerr(1:3)*0)*x(i,5:7)')*kpm-u;
  
% First-Order Filter
 ufilt=phi*ufilt+gam*u_err; 

% Schmidt Trigger   
 if abs(ufilt(1)) < uoff; u(1)=0;end
 if abs(ufilt(2)) < uoff; u(2)=0;end
 if abs(ufilt(3)) < uoff; u(3)=0;end
 if ufilt(1) > uon, u(1)=umax;end; if ufilt(1) < -uon, u(1)=-umax;end
 if ufilt(2) > uon, u(2)=umax;end; if ufilt(2) < -uon, u(2)=-umax;end
 if ufilt(3) > uon, u(3)=umax;end; if ufilt(3) < -uon, u(3)=-umax;end

% Store Torque    
 torque(i,:)=u';

% Integrate 
 f1=dt*euler_quat_reg_thrust_fun(x(i,:),j,invj,u);
 f2=dt*euler_quat_reg_thrust_fun(x(i,:)+0.5*f1',j,invj,u);
 f3=dt*euler_quat_reg_thrust_fun(x(i,:)+0.5*f2',j,invj,u);
 f4=dt*euler_quat_reg_thrust_fun(x(i,:)+f3',j,invj,u);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Error Quaternion
qerr=quat_err(x(:,1:4),[qd(1)*ones(m,1) qd(2)*ones(m,1) qd(3)*ones(m,1) qd(4)*ones(m,1)]);

% Plot Results
clf
subplot(411)
plot(t,qerr(:,1))
axis([0 300 -0.2 0.8])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq1')
subplot(412)
plot(t,qerr(:,2))
axis([0 300 -0.2 0.8])
set(gca,'Ytick',[-0.2 0 0.2 0.4 0.6 0.8])
set(gca,'fontsize',12)
ylabel('dq2')
subplot(413)
plot(t,qerr(:,3))
axis([0 300 -0.2 0.2])
set(gca,'Ytick',[-0.2 -0.1 0 0.1 0.2])
set(gca,'fontsize',12)
ylabel('dq3')
subplot(414)
plot(t,qerr(:,4))
axis([0 300 0 1.2])
set(gca,'Ytick',[0 0.4 0.8 1.2])
set(gca,'fontsize',12)
ylabel('dq4')
xlabel('Time (Sec)')
pause
subplot(311)
plot(t,x(:,5))
set(gca,'fontsize',12)
ylabel('w1 (rad/sec)')
subplot(312)
plot(t,x(:,6))
set(gca,'fontsize',12)
ylabel('w2 (rad/sec)')
subplot(313)
plot(t,x(:,6))
set(gca,'fontsize',12)
ylabel('w3 (rad/sec)')
xlabel('Time (Sex)')
pause
subplot(311)
plot(t,torque(:,1))
axis([0 300 -25 25])
set(gca,'Ytick',[-25 0 25])
set(gca,'fontsize',12)
ylabel('L1 (N-m)')
subplot(312)
plot(t,torque(:,2))
axis([0 300 -25 25])
set(gca,'Ytick',[-25 0 25])
set(gca,'fontsize',12)
ylabel('L2 (N-m)')
subplot(313)
plot(t,torque(:,3))
axis([0 300 -25 25])
set(gca,'Ytick',[-25 0 25])
set(gca,'fontsize',12)
ylabel('L3 (N-m)')
xlabel('Time (Sec)')