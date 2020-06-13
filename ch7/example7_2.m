% This example shows performs a reorientation manuever for the spacecraft
% regulation case involving wheels.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.2

% Written by John L. Crassidis 9/13

% Other Required Routines: wheel_fun.m, quat_err.m

% True Inertia
intrue=[6400 -76.4 -25.6;-76.4 4730 -40;-25.6 -40 8160];

% Time
dt=1;tf=20*60;
t=[0:dt:tf]';
m=length(t);

% Initial Conditions
q0=[sqrt(2)/2;0;0;sqrt(2)/2];
w0=[0.01;0.01;0.01];
wheel0=[0;0;0];

% Pre-Allocate Space
q=zeros(m,4);q(1,:)=q0(:)';
w=zeros(m,3);w(1,:)=w0(:)';
wheel=zeros(m,3);wheel(1,:)=wheel0(:)';
xa=zeros(m,10);xa(1,:)=[q(1,:) w(1,:) wheel(1,:)];

% Desired Quaternion
q_c=[0;0;0;1];
xi_c=[q_c(4) -q_c(3) q_c(2)
         q_c(3)  q_c(4) -q_c(1)
        -q_c(2) q_c(1) q_c(4)
        -q_c(1) -q_c(2) -q_c(3)]; 

% Gains
kp=10;kd=150;

% Main Loop
for i=1:m-1,

dq13=xi_c'*q(i,:)';dq4=q(i,:)*q_c;
wheel_torq=-kp*sign(dq4)*dq13-kd*w(i,:)';

f1=wheel_fun(xa(i,:),intrue,wheel_torq);
f2=wheel_fun(xa(i,:)+0.5*f1'*dt,intrue,wheel_torq);
f3=wheel_fun(xa(i,:)+0.5*f2'*dt,intrue,wheel_torq);
f4=wheel_fun(xa(i,:)+f3'*dt,intrue,wheel_torq);
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion, Angular Velocity and Wheel Speed
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);
wheel(i+1,:)=xa(i+1,8:10);

end

qerr=quat_err(q,kron(q_c',ones(m,1)));

% Plots
clf
subplot(411)
plot(t/60,qerr(:,1));
set(gca,'fontsize',12);
axis([0 20 -0.4 0.8])
set(gca,'ytick',[-0.4 0 0.4 0.8])
ylabel('dq1');
subplot(412)
plot(t/60,qerr(:,2));
set(gca,'fontsize',12);
axis([0 20 -0.04 0.02])
set(gca,'ytick',[-0.04 -0.02 0 0.02])
ylabel('dq2');
subplot(413)
plot(t/60,qerr(:,3));
set(gca,'fontsize',12);
axis([0 20 -0.1 0.2])
set(gca,'ytick',[-0.1 0 0.1 0.2])
ylabel('dq3');
subplot(414)
plot(t/60,qerr(:,4));
set(gca,'fontsize',12);
axis([0 20 0.6 1.2])
set(gca,'ytick',[0.6 0.8 1 1.2])
ylabel('dq4');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,wheel(:,1));
set(gca,'fontsize',12);
axis([0 20 0 240])
set(gca,'ytick',[0 60 120 180 240])
ylabel('h1 (Nms)');
subplot(312)
plot(t/60,wheel(:,2));
set(gca,'fontsize',12);
axis([0 20 -120 40])
set(gca,'ytick',[-120 -80 -40 0 40])
ylabel('h2 (Nms)');
subplot(313)
plot(t/60,wheel(:,3));
set(gca,'fontsize',12);
axis([0 20 0 120])
set(gca,'ytick',[0 40 80 120])
ylabel('h3 (Nms)');
xlabel('Time (Min)');

pause,clf
tot_mom=(intrue*w')'+wheel;
plot(t/60,(tot_mom(:,1).^2+tot_mom(:,2).^2+tot_mom(:,3).^2).^(0.5))
set(gca,'fontsize',12);
ylabel('Total Momentum (N-m-s)');
xlabel('Time (Min)');