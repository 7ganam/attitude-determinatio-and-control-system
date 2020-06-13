% This program simulates the SAMPEX control modes mode
% Written by John L. Crassidis 11/12/12

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis

% Written by John L. Crassidis 11/12

% Other Required Routines: attm.m, crossm.m, eclipse.m, gc2gd.m, julian.m,
% mag_field.m, orbitfun1.m, quat_err.m, solar.m, sun_sensor_sampex_model, wheel_torque_fun 

% Pick Mode to Simulate
mode=input('Desired Mode: 1 = Vertical Pointing, 2 = Orbit Rate Rotation, 3 = Special Pointing ');

% Time Interval
dt=0.5;tf=5*60*60;time=[0:dt:tf]';m=length(time);

% Orbital Elements
mu=398600.4415;rearth=6378.1363;

% Inertia Matrix
in=[10.4 -0.14 -0.73
   -0.14 14.2 -0.36
   -0.73 -0.36 8.9]*1.355817962; % kg-m^2

% Wheel Inertia
iw=3.06e-3*1.355817962; % kg-m^2

% Initial Condition and Other Variables
yr=2011;mth=9;day=16;hr=0;min=0;sec=0;
r0=[3335.973299;2571.763319;-5370.931739];
v0=[3.530941;4.977268;4.566940];
x0=[r0;v0];

% Get True Orbit
x_true=zeros(m,6);
x_true(1,:)=x0(:)';
  for i = 1:m-1
   f1=dt*orbitfun1(x_true(i,:)',mu,rearth);
   f2=dt*orbitfun1(x_true(i,:)'+0.5*f1,mu,rearth);
   f3=dt*orbitfun1(x_true(i,:)'+0.5*f2,mu,rearth);
   f4=dt*orbitfun1(x_true(i,:)'+f3,mu,rearth);
   x_true(i+1,:)=x_true(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
  end
pos=x_true(:,1:3);
vel=x_true(:,4:6);

% Greenwich Mean Sidereal Time from Vallado (third edition), p. 191
jdate=julian(yr,mth,day,hr,min,sec+time);
tdays=jdate-2451545;
jcent=tdays/36525;
gst_sec=67310.54841+(876600*3600+8640184.812866)*jcent+0.093104*jcent.^2-6.2e-6*jcent.^3;
gst=rem(gst_sec,86400)/240;
i_neg=find(gst<0);gst(i_neg)=gst(i_neg)+360;

% Longitude, Latitude and Height of Satellite
[long,lat,height]=gc2gd(pos,yr,mth,day,hr,min,sec,dt,tf);

% Magnetic Field
[mag_i,mag_ecef,mag_ned,ti,bh,dec,dip]=mag_field(pos,yr,mth,day,hr,min,sec,10,dt,1);
mag_i_norm=(mag_i(:,1).^2+mag_i(:,2).^2+mag_i(:,3).^2).^(0.5);

% Sun in Inertial Coordinates and Sun Availability
[sun_i,r_sun]=solar(yr,mth,day,hr,min,sec,dt,tf);
sun_avail=eclipse(pos,[r_sun.*sun_i(:,1) r_sun.*sun_i(:,2) r_sun.*sun_i(:,3)]);

% Angle Between Sun and Magnetic Field Vectors
ang_sun_tam=acos(sun_i(:,1).*mag_i(:,1)./mag_i_norm+sun_i(:,2).*mag_i(:,2)./mag_i_norm+sun_i(:,3).*mag_i(:,3)./mag_i_norm)*180/pi;

% Time to Point y-axis to Sun (wheel is off during this time)
time_sun_point=2; % hours
i_sun_point=time_sun_point*60*60/dt;

% Control Gains for Vertical Pointing and Orbit Rate Rotation Modes
if (mode==1) | (mode==2) 
 zeta=sqrt(2)/2;wn=0.02; 
 kp=in(2,2)*wn^2;
 kd=2*zeta*wn*in(2,2);
 ki=kp*kd/(10*in(2,2));
else
% Control Gain for Special Pointing Mode
 zeta=sqrt(2)/2;wn=0.01; 
 kp=in(2,2)*wn^2;
 kd=2*zeta*wn*in(2,2);
 ki=kp*kd/(10*in(2,2));   
end
int_limit=3e-3/ki;
 
% Integral Control Initial Condition
err_ang_int=0;

% Mag Control Gain
k_mag=7.376e-13;

% Desired Momentum
h_0=0.8135; % N-m-s

% Kalman Gain for Momentum Filter
kal_gain=0.04; % note: use 0.01 if using TRIAD solution

% Velocity Avoidance Minimum Ram Angle
phi_min=80*pi/180;
cos_phi_min=cos(phi_min);

% Initial Conditions for Spacecraft (attitude = I and all momentums = 0)
x=zeros(m,8);x(1,:)=[0 0 0 1 0 0 0 0];
q=zeros(m,4);q(1,:)=x(1,1:4);
h=zeros(m,3);h(1,:)=x(1,5:7);h_est=h(1,:)';
h_wheel=zeros(m,1);h_wheel(1,:)=x(1,8);

% Pre-allocate Space
mag_b=zeros(m,3);sun_b=zeros(m,3);
mag_bm=zeros(m,3);sun_bm=zeros(m,3);
err_ang=zeros(m,1);
wheel_torq=zeros(m,1);
u_mag=zeros(m,3);
w_est=zeros(m,3);
q_est=zeros(m,4);
v_body1=zeros(m,1);

% Measurement Noise Values and Sun Sensor Model Parameters
sig_tam=30; %nT
res_tam=31.25; %nT
sig_sun=0.01; % deg
res_sun=0.5; % deg
t_glass=0.448; % cm
n_refrac=1.4553;
sun_scale=0.002754; % cm/count
sun_bias=-0.350625; % cm

% Main Loop
for i = 1:m-1
    
% Get Attitude Matrix
att_mat=attm(q(i,:));

% Get Body TAM and Sun
mag_b(i,:)=(att_mat*mag_i(i,:)')';
sun_b(i,:)=(att_mat*sun_i(i,:)')';

% Form Measurements
mag_b_noise=mag_b(i,:)+sig_tam*randn(1,3);
mag_bm(i,:)=[round(mag_b_noise(1)/res_tam)*res_tam round(mag_b_noise(2)/res_tam)*res_tam round(mag_b_noise(3)/res_tam)*res_tam];
if sun_avail(i)==1
 sun_bm(i,:)=sun_sensor_sampex_model(sun_b(i,:),n_refrac,t_glass,res_sun,sig_sun)';
else
 sun_bm(i,:)=[0 1 0];
end
 
% q-method Solution
om_mag=[-crossm(mag_bm(i,:)) mag_bm(i,:)';-mag_bm(i,:) 0];
gam_mag=[crossm(mag_i(i,:)) mag_i(i,:)';-mag_i(i,:) 0];
om_sun=[-crossm(sun_bm(i,:)) sun_bm(i,:)';-sun_bm(i,:) 0];
gam_sun=[crossm(sun_i(i,:)) sun_i(i,:)';-sun_i(i,:) 0];   
k_mat=-1/sig_tam^2*om_mag*gam_mag-1/sig_sun^2*om_sun*gam_sun;
[vv,ee]=eig(k_mat);
[ii,jj]=max(diag(ee));
q_est(i,:)=vv(:,jj)'*sign(vv(4,jj));
att_mat_est=attm(q(i,:));

% Velocity in Body Coordinates
v_body=att_mat_est*vel(i,:)'/norm(vel(i,:));
v_body1(i)=v_body(1);

% Vertical Pointing Mode
if mode == 1,
 vel_avoid_flag=1;
 u=att_mat_est*cross(sun_i(i,:)',cross(pos(i,:)',sun_i(i,:)'))/norm(cross(sun_i(i,:)',cross(pos(i,:)',sun_i(i,:)')));
end

% Orbit Rate Rotation Mode
if mode == 2,
 vel_avoid_flag=1;
 orb_normal=cross(pos(i,:),vel(i,:))'/norm(cross(pos(i,:),vel(i,:))');
 w_vec=cross([0;0;1],sun_i(i,:)')/norm(cross([0;0;1],sun_i(i,:)'));
 an_vec=cross([0;0;1],orb_normal)/norm(cross([0;0;1],orb_normal));
 nmp_vec=cross(orb_normal,an_vec);

 sin_alpha=-pos(i,:)*an_vec/norm(pos(i,:));
 cos_alpha=pos(i,:)*nmp_vec/norm(pos(i,:));

 if sun_i(i,:)*orb_normal<0
  u_i=cos_alpha*cross(sun_i(i,:)',w_vec)+sin_alpha*w_vec;
 else
  u_i=cos_alpha*cross(sun_i(i,:)',w_vec)-sin_alpha*w_vec; 
 end
 
 u=att_mat_est*u_i;
end

% Special Pointing Mode    
if mode == 3,
 if (mag_i_norm(i)<=3e4) & (vel(i,:)*cross(sun_i(i,:)',mag_i(i,:)')<=0)
  u=1/sqrt(mag_bm(i,1)^2+mag_bm(i,3)^2)*[mag_bm(i,3);0;-mag_bm(i,1)];
  vel_avoid_flag=0;
 end
 if (mag_i_norm(i)<=3e4) & (vel(i,:)*cross(sun_i(i,:)',mag_i(i,:)')>0)
  u=-1/sqrt(mag_bm(i,1)^2+mag_bm(i,3)^2)*[mag_bm(i,3);0;-mag_bm(i,1)];
  vel_avoid_flag=0;
 end
 if (mag_i_norm(i)>3e4) & (lat(i)<=0)
  u=1/sqrt(mag_bm(i,1)^2+mag_bm(i,3)^2)*[mag_bm(i,1);0;mag_bm(i,3)];
  vel_avoid_flag=1;
 end
 if (mag_i_norm(i)>3e4) & (lat(i)>0)
  u=-1/sqrt(mag_bm(i,1)^2+mag_bm(i,3)^2)*[mag_bm(i,1);0;mag_bm(i,3)];
  vel_avoid_flag=1;
 end
 
end

% Velocity Avoidance
w_f=cross(u,sun_bm(i,:)')/norm(cross(u,sun_bm(i,:)'));
if (v_body'*u <= cos_phi_min) | (vel_avoid_flag == 0)
 sin_theta_va=0;cos_theta_va=1;
else
 v_f=[sun_bm(i,:);w_f';u']*v_body;
 term1=v_f(2)*cos_phi_min;
 arg=v_f(2)^2+v_f(3)^2-cos_phi_min^2;
 if (arg < 0)
  arg=0;
 end
 term2=abs(v_f(3))*sqrt(arg);
 num1=term1+term2;
 num2=term1-term2;
 den=v_f(2)^2+v_f(3)^2;
 if (v_f(2) >= 0)
  sin_theta_va=num2/den;
 else
  sin_theta_va=num1/den;
 end
 if (sin_theta_va > 1)
  sin_theta_va=1;
 elseif (sin_theta_va < -1)
  sin_theta_va=-1;
 end
 theta=asin(sin_theta_va); 
 sin_theta_va=sin(theta);cos_theta_va=cos(theta);
end
u_ram=sin_theta_va*w_f+cos_theta_va*u;

% Error Angle (note: change in angle for low-field and large maneuver)
if (vel_avoid_flag == 0) & (abs(atan2(-u_ram(1),u_ram(3))) > 2.5);
  err_ang(i)=atan2(-sign(v_body(1))*u_ram(1),sign(v_body(1))*u_ram(3));
else
  err_ang(i)=atan2(-u_ram(1),u_ram(3));  
end 
 
% Control Law
if i > 1
xi=[q_est(i-1,4)*eye(3)+crossm(q_est(i-1,1:3)');-q_est(i-1,1:3)];
w_est(i,:)=2*(xi'*q_est(i,:)')'/dt;

% Wheel Torque
if (i < i_sun_point) % keep wheel off until y-axis is pointed towards Sun
 wheel_torq(i)=0; % note: h_wheel(1) = 0 so h_wheel remains at zero
 tam_flag=1;
elseif (ang_sun_tam(i) < 5) & (sun_avail(i)==1)
 wheel_torq(i)=wheel_torq(i-1);
 tam_flag=1;
elseif (ang_sun_tam(i) < 40) & (sun_avail(i)==0) 
 wheel_torq(i)=wheel_torq(i-1);
 tam_flag=0; 
else
 err_ang_int=err_ang_int+err_ang(i)*dt;
 if (abs(err_ang_int) > int_limit) 
  err_ang_int=int_limit*sign(err_ang_int);
 end
 err_ang_rate=(err_ang(i)-err_ang(i-1))/dt;
 wheel_torq(i)=(kp*err_ang(i)+kd*err_ang_rate+ki*err_ang_int);
 tam_flag=1;
end
% Mag Torque
h_derived=in*w_est(i,:)'+h_wheel(i)*[0;1;0];
h_pred=-crossm(w_est(i-1,:))*h_est*dt+h_est+u_mag(i-1,:)'*dt;
h_est=(1-kal_gain)*h_pred+kal_gain*h_derived;
delta_h=2*h_est-h_0*([0;1;0]+sun_bm(i,:)');
m_mag=k_mag*cross(delta_h,mag_bm(i,:)');
u_mag(i,:)=cross(m_mag',mag_bm(i,:))*tam_flag;
end

% Integration Routine
f1=wheel_torque_fun(x(i,:),in,wheel_torq(i),u_mag(i,:));
f2=wheel_torque_fun(x(i,:)+0.5*f1'*dt,in,wheel_torq(i),u_mag(i,:));
f3=wheel_torque_fun(x(i,:)+0.5*f2'*dt,in,wheel_torq(i),u_mag(i,:));
f4=wheel_torque_fun(x(i,:)+f3'*dt,in,wheel_torq(i),u_mag(i,:));
x(i+1,:)=x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion, Angular Velocity and Wheel Speed
q(i+1,:)=x(i+1,1:4);
h(i+1,:)=x(i+1,5:7);
h_wheel(i+1)=x(i+1,8);
    
end

% Set Last Point to Previous Point (just for plots)
sun_b(i+1,:)=sun_b(i,:);
err_ang(i+1)=err_ang(i);
q_est(i+1,:)=q_est(i,:);
u_mag(i+1,:)=u_mag(i,:);

% True Angular Velocity
w_true=(inv(in)*(h'-[zeros(1,m);h_wheel';zeros(1,m)]))';

% Estimation Angle Errors
q_est=q_est.*kron([1 1 1 1],sign(q_est(:,4).*q(:,4)));
qerr=quat_err(q_est,q);erre=qerr(:,1:3)*2*180/pi;

% Plot Results
clf
load topo
imageHndl=imagesc(topo);
set(gca,'Ydir','normal','XLim',[0 360],'Ylim',[0 180])
colormap((topomap1+white)/2);
hold on;
set(gca,'Xtick',[0:30:360]')
set(gca,'XtickLabel',[0 30 60 90 120 150 180 -150 -120 -90 -60 -30 0])
set(gca,'Ytick',[0:10:180]')
set(gca,'YtickLabel',[-90,-80,-70,-60,-50,-40,-30,-20,-10,0, ...
10 20 30 40 50 60 70 80 90])
grid
xlabel('Longitude (Deg)')
ylabel('Latitude (Deg)')
plot(rem(unwrap(rem(long,360)+360),360),lat+90,'r.')
hold off;
clear topo

pause
% Sun Angle
plot(time/60/60,acos(sun_b(:,2))*180/pi)
set(gca,'fontsize',12)
axis([0 tf/60/60 -20 100])
ylabel('Sun Angle Error (Deg)')
xlabel('Time (Hr)')

pause
% Error Angle
plot(time/60/60,err_ang*180/pi)
set(gca,'fontsize',12)
ylabel('Angle Error (Deg)')
xlabel('Time (Hr)')

pause
% Wheel Momentum
plot(time/60/60,h_wheel)
set(gca,'fontsize',12)
axis([0 tf/60/60 -0.2 1.2])
set(gca,'ytick',[-0.2 0 0.2 0.4 0.6 0.8 1 1.2])
ylabel('Wheel Momentum (Nms)')
xlabel('Time (Hr)')

pause
% Magnetic Torques
subplot(311)
plot(time/60/60,u_mag(:,1))
axis([0 tf/60/60 -1e-3 1e-3])
set(gca,'ytick',[-1e-3 -0.5e-3 0 0.5e-3 1e-3])
set(gca,'fontsize',12)
ylabel('x (Nm)')
subplot(312)
plot(time/60/60,u_mag(:,2))
axis([0 tf/60/60 -1e-3 1e-3])
set(gca,'ytick',[-1e-3 -0.5e-3 0 0.5e-3 1e-3])
set(gca,'fontsize',12)
ylabel('y (Nm)')
subplot(313)
plot(time/60/60,u_mag(:,3))
axis([0 tf/60/60 -1e-3 1e-3])
set(gca,'ytick',[-1e-3 -0.5e-3 0 0.5e-3 1e-3])
set(gca,'fontsize',12)
ylabel('z (Nm)')
xlabel('Time (Hr)')