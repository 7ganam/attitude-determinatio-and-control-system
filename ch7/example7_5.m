% This example shows performs detumbles a spacecraft using
% magnetic torquers

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.5

% Written by John L. Crassidis 9/13

% Other Required Routines: att.m, euler_quat_fun.m, julian.m, mag_field.m, mytle.m, orbitfun.m, tle.txt

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

% Time and Pre-allocation of Variables
dt=5;tf=3*60*60;
t=[0:dt:tf]';m=length(t);
x_pos_vel=zeros(m,6);x_pos_vel(1,:)=[r0' v0'];
x=zeros(m,7);u=zeros(m,3);command_dipole=zeros(m,3);
q0=[sqrt(2)/2;0;0;sqrt(2)/2];
w0=[0.01;0.01;0.01];
x(1,:)=[q0' w0'];

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

% Gain for Controller
gain=2*(2*pi/orb_per)*(1+sin(dip))*min_prin_mom;

% Main Loop
for i=1:m-1
    
 b_body=attm(x(i,1:4)')*b_eci(i,:)';b_body_n=b_body/norm(b_body);
 u(i,:)=(-gain(i)*(eye(3)-b_body_n*b_body_n')*x(i,5:7)')';
 command_dipole(i,:)=(-gain(i)/norm(b_body)*cross(b_body_n,(eye(3)-b_body_n*b_body_n')*x(i,5:7)'))'*1e9;% b_body is in nT. Need to convert to T.;
 
 f1=dt*euler_quat_fun(x(i,:),in,u(i,:)');
 f2=dt*euler_quat_fun(x(i,:)+0.5*f1',in,u(i,:)');
 f3=dt*euler_quat_fun(x(i,:)+0.5*f2',in,u(i,:)');
 f4=dt*euler_quat_fun(x(i,:)+f3',in,u(i,:)');
 x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

b_body=attm(x(i+1,1:4)')*b_eci(i+1,:)';b_body_n=b_body/norm(b_body);
u(i+1,:)=(-gain(i+1)*(eye(3)-b_body_n*b_body_n')*x(i+1,5:7)')';
command_dipole(i+1,:)=(-gain(i+1)/norm(b_body)*cross(b_body_n,(eye(3)-b_body_n*b_body_n')*x(i+1,5:7)'))';

clf
subplot(311)
plot(t/60,x(:,5)*180/pi);
set(gca,'fontsize',12);
axis([0 180 -0.6 0.6])
set(gca,'ytick',[-0.6 -0.3 0 0.3 0.6])
ylabel('w1 (Deg/Sec)');
subplot(312)
plot(t/60,x(:,6)*180/pi);
set(gca,'fontsize',12);
axis([0 180 -0.3 0.9])
set(gca,'ytick',[-0.3 0 0.3 0.6 0.9])
ylabel('w2 (Deg/Sec)');
subplot(313)
plot(t/60,x(:,7)*180/pi);
set(gca,'fontsize',12);
axis([0 180 -0.3 0.6])
set(gca,'ytick',[-0.6 -0.3 0 0.3 0.6])
ylabel('w3 (Deg/Sec)');
xlabel('Time (Min)');

pause

subplot(311)
plot(t/60,u(:,1));
set(gca,'fontsize',12);
axis([0 180 -0.08 0.08])
set(gca,'ytick',[-0.08 -0.04 0 0.04 0.08])
ylabel('L1 (N-m)');
subplot(312)
plot(t/60,u(:,2));
set(gca,'fontsize',12);
axis([0 180 -0.2 0.1])
set(gca,'ytick',[-0.2 -0.1 0 0.1])
ylabel('L2 (N-m)');
subplot(313)
plot(t/60,u(:,3));
set(gca,'fontsize',12);
axis([0 180 -0.12 0.04])
set(gca,'ytick',[-0.12 -0.08 -0.04 0 0.04])
ylabel('L3 (N-m)');
xlabel('Time (Min)');
