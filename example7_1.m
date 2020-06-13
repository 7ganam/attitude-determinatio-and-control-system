% This example shows performs a reorientation manuever for the spacecraft
% regulation case involving external torques.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.1

% Written by John L. Crassidis 9/13

% Other Required Routines: euler_quat_reg_fun.m, quat_err.m

% Time and Inertia Matrix
dt=1;tf=300;
t=[0:dt:tf]';m=length(t);
j=diag([10000 9000 12000]);invj=inv(j);

% Initial and Desired Quaternion, and Initial Rate
qd=[0;0;0;1];
q0=[0.685;0.695;0.153;0.153];q0=q0/norm(q0);
w0=[0.53;0.53;0.053]*pi/180;
x=zeros(m,7);x(1,:)=[q0' w0'];

% Gains
kp=50;kd=500;

% Main Loop
for i = 1:m-1
 f1=dt*euler_quat_reg_fun(x(i,:),j,invj,qd,kp,kd);
 f2=dt*euler_quat_reg_fun(x(i,:)+0.5*f1',j,invj,qd,kp,kd);
 f3=dt*euler_quat_reg_fun(x(i,:)+0.5*f2',j,invj,qd,kp,kd);
 f4=dt*euler_quat_reg_fun(x(i,:)+f3',j,invj,qd,kp,kd);
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Error Quaternion
qerr=quat_err(x(:,1:4),[qd(1)*ones(m,1) qd(2)*ones(m,1) qd(3)*ones(m,1) qd(4)*ones(m,1)]);

% Torque
%torque=-kp*qerr(:,1:3)-kd*x(:,5:7);
torque=-kp*kron(sign(qerr(:,4)),[1 1 1]).*qerr(:,1:3)-kd*x(:,5:7);

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
xlabel('Time (Sec)')
pause
subplot(311)
plot(t,torque(:,1))
set(gca,'fontsize',12)
ylabel('L1 (N-m)')
subplot(312)
plot(t,torque(:,2))
set(gca,'fontsize',12)
ylabel('L2 (N-m)')
subplot(313)
plot(t,torque(:,3))
set(gca,'fontsize',12)
ylabel('L3 (N-m)')
xlabel('Time (Sec)')
