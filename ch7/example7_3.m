% This program tracks an attitude trajectory using sliding mode control.
% The WMAP spacecraft is the model used for the simulation.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.3

% Other Required Routines: wheel_track_fun.m, sat1.m

% Written by John L. Crassidis 8/13

% True Inertia 
intrue=[399,-2.81,-1.31;-2.81,377,2.54;-1.31,2.54,377];

% Assumed Inertia
in=[380,-2.81,-1.31;-2.81,360,2.54;-1.31,2.54,340];

% Time
dt=0.1;tf=3600;
t=[0:dt:tf]';
m=length(t);

% Pre-Allocate Space
q=zeros(m,4);
q_d=zeros(m,4);
w_d=zeros(m,3);
w_dd=zeros(m,3);
w=zeros(m,3);
u=zeros(m,3);
wheel=zeros(m,3);
slide=zeros(m,3);
xa=zeros(m,10);
i1000=0;

% Disturbance
dist=[0.005*sin(0.05*t) 0.003*ones(m,1) 0.005*cos(0.05*t)]; 

% Set null=1 for Tracking and null=0 for Regulation 
null=1;

% Desired Euler Quantities
phir_d=1*(2*pi/3600)*null;
theta_d=22.5*pi/180*null;
psir_d=0.464*2*pi/60*null;

% Get Euler Angles (zero conditions)
phi_d=phir_d.*t*null;
psi_d=psir_d.*t*null;

% Get Initial Desired Quaternion and Angular Velocity
q_d(1,1)=sin(theta_d/2).*cos((phi_d(1)-psi_d(1))/2);
q_d(1,2)=sin(theta_d/2).*sin((phi_d(1)-psi_d(1))/2);
q_d(1,3)=cos(theta_d/2).*sin((phi_d(1)+psi_d(1))/2);
q_d(1,4)=cos(theta_d/2).*cos((phi_d(1)+psi_d(1))/2);

w_d(1,1)=sin(theta_d).*sin(psi_d(1)).*phir_d;
w_d(1,2)=sin(theta_d).*cos(psi_d(1)).*phir_d;
w_d(1,3)=cos(theta_d).*phir_d+psir_d;


% Initial Quaternion and Rate
qc_d=[0 -q_d(1,3) q_d(1,2)
   q_d(1,3) 0 -q_d(1,1)
   -q_d(1,2) q_d(1,1) 0];
xiq_d=[q_d(1,4)*eye(3)+qc_d;-q_d(1,1:3)];
q(1,:)=([xiq_d q_d(1,:)']*[0;0;sin(60/2*pi/180);cos(60/2*pi/180)])';
w(1,:)=w_d(1,:)*0;
xa(1,:)=[q(1,:) w(1,:) 0 0 0];

% Gains
k=0.015;g=0.15*eye(3);eps=0.01;

% Get Torque, Sliding Manifold, Desired Quaternion and Angular Velocity at Initial Time
[ff,torq,ss,qq,ww]=wheel_track_fun(xa(1,:),t(1),theta_d,phir_d,psir_d,dist(1,:),in,intrue,k,g,eps);
u(1,:)=torq(:)';
slide(1,:)=ss(:)';
q_d(1,:)=qq(:)';
w_d(1,:)=ww(:)';

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

f1=wheel_track_fun(xa(i,:),t(i),theta_d,phir_d,psir_d,dist(i,:),in,intrue,k,g,eps);
f2=wheel_track_fun(xa(i,:)+0.5*f1'*dt,t(i)+0.5*dt,theta_d,phir_d,psir_d,dist(i,:),in,intrue,k,g,eps);
f3=wheel_track_fun(xa(i,:)+0.5*f2'*dt,t(i)+0.5*dt,theta_d,phir_d,psir_d,dist(i,:),in,intrue,k,g,eps);
f4=wheel_track_fun(xa(i,:)+f3'*dt,t(i)+dt,theta_d,phir_d,psir_d,dist(i,:),in,intrue,k,g,eps);
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Get Torque, Sliding Manifold, Desired Quaternion and Angular Velocity
[ff,torq,ss,qq,ww]=wheel_track_fun(xa(i+1,:),t(i+1,:),theta_d,phir_d,psir_d,dist(i+1,:),in,intrue,k,g,eps);
u(i+1,:)=torq(:)';
slide(i+1,:)=ss(:)';
q_d(i+1,:)=qq(:)';
w_d(i+1,:)=ww(:)'; 

% Actual Quaterion, Angular Velocity and Wheel Speed
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);
wheel(i+1,:)=xa(i+1,8:10);

end

% Get Roll, Pitch and Yaw (3-1-3 sequence)
phi=atan2(2*(q(:,1).*q(:,3)+q(:,2).*q(:,4)),-2*(q(:,2).*q(:,3)-q(:,1).*q(:,4)));
theta=acos(-q(:,1).^2-q(:,2).^2+q(:,3).^2+q(:,4).^2);
psi=atan2(2*(q(:,1).*q(:,3)-q(:,2).*q(:,4)),2*(q(:,2).*q(:,3)+q(:,1).*q(:,4)));
werr=(w_d-w)*180/pi;
erre=[phi_d-unwrap(phi) theta_d-theta psi_d-(unwrap(psi)+2*pi)]*180/pi;
erre(:,3)=erre(:,3)+360;

% Plots
clf
subplot(311)
plot(t/60,erre(:,1));
set(gca,'fontsize',12);
axis([0 60 -2 2])
set(gca,'ytick',[-2 -1 0 1 2])
ylabel('Roll (Deg)');
subplot(312)
plot(t/60,erre(:,2));
set(gca,'fontsize',12);
axis([0 60 -4 4])
set(gca,'ytick',[-4 -2 0 2 4])
ylabel('Pitch (Deg)');
subplot(313)
plot(t/60,erre(:,3));
set(gca,'fontsize',12);
axis([0 60 -60 30])
set(gca,'ytick',[-60 -30 0 30])
ylabel('Yaw (Deg)');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,werr(:,1));
set(gca,'fontsize',12);
axis([0 60 -0.02 0.02])
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02])
ylabel('dw1 (Deg/Sec)');
subplot(312)
plot(t/60,werr(:,2));
set(gca,'fontsize',12);
axis([0 60 -0.02 0.02])
set(gca,'ytick',[-0.02 -0.01 0 0.01 0.02])
ylabel('dw2 (Deg/Sec)');
subplot(313)
plot(t/60,werr(:,3));
set(gca,'fontsize',12);
axis([0 60 -0.2 0.4])
set(gca,'ytick',[-0.2 0 0.2 0.4])
ylabel('dw3 (Deg/Sec)');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,wheel(:,1));
set(gca,'fontsize',12);
axis([0 20 -0.4 0.4])
set(gca,'ytick',[-0.4 -0.2 0 0.2 0.4])
ylabel('h1 (Nms)');
subplot(312)
plot(t/60,wheel(:,2));
set(gca,'fontsize',12);
axis([0 20 -0.4 0.2])
set(gca,'ytick',[-0.4 -0.2 0 0.2])
ylabel('h2 (Nms)');
subplot(313)
plot(t/60,wheel(:,3));
set(gca,'fontsize',12);
axis([0 20 -20 -16])
set(gca,'ytick',[-20 -19 -18 -17 -16])
ylabel('h3 (Nms)');
xlabel('Time (Min)');

pause

d_norm=(dist(:,1).^2+dist(:,2).^2+dist(:,3).^2).^(0.5);
d_norm_max=max(d_norm);
bound=norm(eps*inv(in*g),'fro')*d_norm_max*ones(m,1);
slide_norm=(slide(:,1).^2+slide(:,2).^2+slide(:,3).^2).^(0.5);
clf
plot(t/60,slide_norm,t/60,bound)
set(gca,'fontsize',12);
axis([0 60 0 3e-6])
ylabel('Slide Norm and Bound');
xlabel('Time (Min)');

pause

% Scan Plot
e1=cos(theta_d).*cos(psi_d).*cos(phi_d)-cos(theta_d).^2.*sin(psi_d).*sin(phi_d)+sin(theta_d).^2.*sin(phi_d);
e2=cos(theta_d).*cos(psi_d).*sin(phi_d)+cos(theta_d).^2.*sin(psi_d).*cos(phi_d)-sin(theta_d).^2.*cos(phi_d);
e3=cos(theta_d).*sin(theta_d).*(sin(psi_d)+1);
x=e1./(1+e3);y=e2./(1+e3);
plot(x,y);
set(gca,'fontsize',12);
set(gca,'XTicklabels',[])
set(gca,'YTicklabels',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'DataAspectRatio',[1 1 1])
axis([-1.05 1.05 -1.05 1.05])
