% This program regulates a spacecraft using wheels and dumps
% momentum by about half in about 3 hours.
% The TRMM spacecraft is the model used for the simulation.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.6

% Written by John L. Crassidis 9/13

% Other Required Routines: attm.m, julian.m, mag_field.m, mytle.m, quat_error.m, orbitfun.m, wheel_torque_fun.m, tle.txt

% Time 
dt=5;tf=3*60*60;
t=[0:dt:tf]';m=length(t);

% Mu value
mu=398600.64;

% Inertia
in=[6400 -76.4 -25.6;-76.4 4730 -40;-25.6 -40 8160];
[vv,prin_mom]=eig(in);min_prin_mom=min(diag(prin_mom));

% TRMM Spacecraft (Use ecc=0)
tlef='tle.txt';tsat ='TRMM';
sat=mytle(tlef,tsat,mu)
yr=sat.epoch(1);mth=sat.epoch(2);day=sat.epoch(3);hr=sat.epoch(4);minute=sat.epoch(5);sec=sat.epoch(6);
a=sat.oe(1);ecc=sat.oe(2)*0;inc=sat.oe(3)*pi/180;big_omega=sat.oe(4)*pi/180;w=sat.oe(5)*pi/180;big_m=sat.oe(6)*pi/180;
b_star=sat.bStar;b_coef=1/(12.741621*b_star);

% Solve Kepler's Equation
% Initial Guess for Eccentric Anomaly
big_e=big_m;
delta_e=10;eps=1e-15;
max_iter=100;count=1;

while abs(delta_e) > eps
  delta_e=(big_m-(big_e-ecc*sin(big_e)))/(1-ecc*cos(big_e));
  big_e=big_e+delta_e;
  count = count + 1;
  if count == max_iter, break, disp(' Maximum Number of Iterations Achieved'), end
end

% Get Initial Position and Velocity
rmag=a*(1-ecc*cos(big_e));
rp=a*(cos(big_e)-ecc);rq=a*sqrt(1-ecc^2)*sin(big_e);
vp=-sqrt(mu*a)/rmag*sin(big_e);vq=sqrt(mu*a*(1-ecc^2))/rmag*cos(big_e);
c11=cos(big_omega)*cos(w)-sin(big_omega)*sin(w)*cos(inc);
c12=-cos(big_omega)*sin(w)-sin(big_omega)*cos(w)*cos(inc);
c21=sin(big_omega)*cos(w)+cos(big_omega)*sin(w)*cos(inc);
c22=-sin(big_omega)*sin(w)+cos(big_omega)*cos(w)*cos(inc);
c31=sin(w)*sin(inc);
c32=cos(w)*sin(inc);
r1=c11*rp+c12*rq;r2=c21*rp+c22*rq;r3=c31*rp+c32*rq;
v1=c11*vp+c12*vq;v2=c21*vp+c22*vq;v3=c31*vp+c32*vq;
r0=[r1;r2;r3];
v0=[v1;v2;v3];

x_pos_vel=zeros(m,6);x_pos_vel(1,:)=[r0' v0'];

% Orbit Period
orb_per=2*pi/sqrt(mu)*(norm(r0)^(3/2));

% Get Orbit
for i=1:m-1     
 f1=dt*orbitfun(x_pos_vel(i,:),mu);
 f2=dt*orbitfun(x_pos_vel(i,:)+0.5*f1',mu);
 f3=dt*orbitfun(x_pos_vel(i,:)+0.5*f2',mu);
 f4=dt*orbitfun(x_pos_vel(i,:)+f3',mu);
 x_pos_vel(i+1,:)=x_pos_vel(i,:)+1/6*(f1'+2*f2'+2*f3'+f4'); 
end

% Magnetic Field Reference
[b_eci,b_ecef,b_ned,ti,bh,dec,dip]=mag_field(x_pos_vel(:,1:3),yr,mth,day,hr,minute,sec,10,dt,1);
dip=dip*pi/180;

% Initial Quaternion and Angular Rate Conditions
q0=[sqrt(2)/2;0;0;sqrt(2)/2];
w0=[0.01;0.01;0.01];
wheel0=[0;0;0];

% Desired Angular Velocity and Quaternion
% Note that a nonzero w_c can be used if desired
w_c=[zeros(m,1) -2*pi/orb_per*ones(m,1)*0 zeros(m,1)];
q_c=zeros(m,4);q_c(1,:)=[0 0 0 1]*0;q_c(:,4)=ones(m,1);
% for i = 1:m-1
%  wc_norm=norm(w_c(i,:));
%  co=cos(0.5*wc_norm*dt);
%  si=sin(0.5*wc_norm*dt);
%  n1=w_c(i,1)/wc_norm;n2=w_c(i,2)/wc_norm;n3=w_c(i,3)/wc_norm;
%  qw1=n1*si;qw2=n2*si;qw3=n3*si;qw4=co;
%  om=[qw4  qw3 -qw2 qw1;-qw3  qw4  qw1 qw2;qw2 -qw1  qw4 qw3;-qw1 -qw2 -qw3 qw4];
%  q_c(i+1,1:4)=(om*q_c(i,1:4)')';
% end

% Pre-Allocate Space
q=zeros(m,4);q(1,:)=q0(:)';
w=zeros(m,3);w(1,:)=w0(:)';
wheel=zeros(m,3);wheel(1,:)=wheel0(:)';
u=zeros(m,3);command_dipole=zeros(m,3);
xa=zeros(m,10);xa(1,:)=[q(1,:) w(1,:) wheel(1,:)];

% Gains for Wheels
kp=10;kd=150;

% Gain for Magnetic Controller
gain=0.0001;

% Main Loop
for i=1:m-1,

% Wheel Controller
xi_c=[q_c(i,4) -q_c(i,3) q_c(i,2)
         q_c(i,3)  q_c(i,4) -q_c(i,1)
        -q_c(i,2) q_c(i,1) q_c(i,4)
        -q_c(i,1) -q_c(i,2) -q_c(i,3)]; 
dq13=xi_c'*q(i,:)';dq4=q(i,:)*q_c(i,:)';
wheel_torq=-kp*sign(dq4)*dq13-kd*(w(i,:)-w_c(i,:))';

% Magnetic Controller
b_body=attm(xa(i,1:4)')*b_eci(i,:)';b_body_n=b_body/norm(b_body);
if i < 360
 u(i,:)=[0;0;0];command_dipole(i,:)=[0;0;0];
else
 u(i,:)=(-gain*(eye(3)-b_body_n*b_body_n')*xa(i,8:10)')';
 command_dipole(i,:)=(-gain/norm(b_body)*cross(b_body_n,xa(i,8:10)'))'*1e9;% b_body is in nT. Need to convert to T.;
end

f1=wheel_torque_fun(xa(i,:),in,wheel_torq,u(i,:));
f2=wheel_torque_fun(xa(i,:)+0.5*f1'*dt,in,wheel_torq,u(i,:));
f3=wheel_torque_fun(xa(i,:)+0.5*f2'*dt,in,wheel_torq,u(i,:));
f4=wheel_torque_fun(xa(i,:)+f3'*dt,in,wheel_torq,u(i,:));
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion, Angular Velocity and Wheel Speed
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);
wheel(i+1,:)=xa(i+1,8:10);

end

qerr=quat_err(q,q_c);

% Plots
clf
subplot(411)
plot(t/60,qerr(:,1));
set(gca,'fontsize',12);
axis([0 180 -0.4 0.8])
set(gca,'ytick',[-0.4 0 0.4 0.8])
ylabel('dq1');
subplot(412)
plot(t/60,qerr(:,2));
set(gca,'fontsize',12);
axis([0 180 -0.04 0.02])
set(gca,'ytick',[-0.04 -0.02 0 0.02])
ylabel('dq2');
subplot(413)
plot(t/60,qerr(:,3));
set(gca,'fontsize',12);
axis([0 180 -0.1 0.2])
set(gca,'ytick',[-0.1 0 0.1 0.2])
ylabel('dq3');
subplot(414)
plot(t/60,qerr(:,4));
set(gca,'fontsize',12);
axis([0 180 0.6 1.2])
set(gca,'ytick',[0.6 0.8 1 1.2])
ylabel('dq4');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,wheel(:,1));
set(gca,'fontsize',12);
axis([0 180 -80 240])
set(gca,'ytick',[-80 0 80 160 240])
ylabel('h1 (Nms)');
subplot(312)
plot(t/60,wheel(:,2));
set(gca,'fontsize',12);
axis([0 180 -120 40])
set(gca,'ytick',[-120 -80 -40 0 40])
ylabel('h2 (Nms)');
subplot(313)
plot(t/60,wheel(:,3));
set(gca,'fontsize',12);
axis([0 180 -80 160])
set(gca,'ytick',[-80 0 80 160])
ylabel('h3 (Nms)');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,command_dipole(:,1));
set(gca,'fontsize',12);
axis([0 180 -4e2 4e2])
set(gca,'ytick',[-4e2 -2e2 0 2e2 4e2])
ylabel('m1 (A-m2)','Interpreter','none');
subplot(312)
plot(t/60,command_dipole(:,2));
set(gca,'fontsize',12);
axis([0 180 -4e2 2e2])
set(gca,'ytick',[-4e2 -2e2 0 2e2])
ylabel('m2 (A-m2)','Interpreter','none');
subplot(313)
plot(t/60,command_dipole(:,3));
set(gca,'fontsize',12);
axis([0 180 -4e2 4e2])
set(gca,'ytick',[-4e2 -2e2 0 2e2 4e2])
ylabel('m3 (A-m2)','Interpreter','none');
xlabel('Time (Min)','Interpreter','none');

pause,clf
tot_mom=(in*w')'+wheel;
plot(t/60,(tot_mom(:,1).^2+tot_mom(:,2).^2+tot_mom(:,3).^2).^(0.5))
set(gca,'fontsize',12);
axis([0 180 60 120])
ylabel('Total Momentum (Nms)');xlabel('Time (Min)');

