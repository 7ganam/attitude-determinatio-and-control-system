% This program simulates the WMAP tracker using the conroller from the 
% following paper: Mayhew, C.G., Sanfelice, R.G., and Teel, A.R., 
% "Quaternion-Based Hybrid Control for Robust Global Attitude Tracking," 
% IEEE Transactions on Automatic Control, Vol. 56, No. 11, Nov. 2011,
% pp. 2555-2566.

% Fundamentals of Spacecraft Attitude Determination and Control by Markley and Crassidis
% Example 7.8

% Written by John L. Crassidis 9/13

% Other Required Routines: euler_wmap_fun.m

% True Inertia 
intrue=[399,-2.81,-1.31;-2.81,377,2.54;-1.31,2.54,377];

% Assumed Inertia (same as true one)
in=intrue;

% Time
dt=0.1;tf=3600;
t=[0:dt:tf]';
m=length(t);

% Pre-Allocate Space
q=zeros(m,4);
qm=zeros(m,4);
q_d=zeros(m,4);
w_d=zeros(m,3);
w_dd=zeros(m,3);
w=zeros(m,3);
wm=zeros(m,3);
u=zeros(m,3);
h_store=zeros(m,1);
torq_norm=zeros(m,1);
werr_hybrid=zeros(m,3);
werr_dis=zeros(m,3);
werr_unwind=zeros(m,3);
xa=zeros(m,7);
i1000=0;

% Disturbance (if desired)
dist=[0.005*sin(0.05*t) 0.003*ones(m,1) 0.005*cos(0.05*t)]*0; 

% Set null=1 for Tracking and null=0 for Regulation 
null=1;

% Desired Euler Quantities
phir_d=1*(2*pi/3600)*null;
theta_d=22.5*pi/180*null;
psir_d=0.464*2*pi/60*null;

% Get Euler Angles (zero conditions)
phi_d=phir_d.*t*null;
psi_d=psir_d.*t*null;

% Get Desired Quaternion and Angular Velocity
q_d(:,1)=sin(theta_d/2).*cos((phi_d-psi_d)/2);
q_d(:,2)=sin(theta_d/2).*sin((phi_d-psi_d)/2);
q_d(:,3)=cos(theta_d/2).*sin((phi_d+psi_d)/2);
q_d(:,4)=cos(theta_d/2).*cos((phi_d+psi_d)/2);

w_d(:,1)=sin(theta_d).*sin(psi_d(1)).*phir_d;
w_d(:,2)=sin(theta_d).*cos(psi_d(1)).*phir_d;
w_d(:,3)=cos(theta_d).*phir_d+psir_d;

w_dd(:,1)=psir_d.*sin(theta_d).*cos(psi_d).*phir_d;
w_dd(:,2)=-psir_d.*sin(theta_d).*sin(psi_d).*phir_d;
w_dd(:,3)=zeros(m,1);

% Initial Quaternion and Rate
qc_d=[0 -q_d(1,3) q_d(1,2)
   q_d(1,3) 0 -q_d(1,1)
   -q_d(1,2) q_d(1,1) 0];
xiq_d=[q_d(1,4)*eye(3)+qc_d;-q_d(1,1:3)];
q(1,:)=([xiq_d q_d(1,:)']*[0;0;sin(180/2*pi/180);cos(180/2*pi/180)])';
w(1,:)=-w_d(1,:)*5;
xa(1,:)=[q(1,:) w(1,:)];

% Gains for Controller
kp=5;kd=3;

% Get Torque
% Measurements
sigv=sqrt(10)*1e-7;
w_noise=sqrt(dt)*sigv*randn(m,3);
sig_star=0.5e-3*pi/180/3;
qv=sig_star*randn(m,3);
q_noise_un=[qv ones(m,1)];
q_noise_norm=(q_noise_un(:,1).^2+q_noise_un(:,2).^2+q_noise_un(:,3).^2+q_noise_un(:,4).^2).^(0.5);
q_noise=q_noise_un./[q_noise_norm q_noise_norm q_noise_norm q_noise_norm];
qc=[0 -q(1,3) q(1,2)
   q(1,3) 0 -q(1,1)
   -q(1,2) q(1,1) 0];
wm(1,:)=w(1,:)+w_noise(1,:);
qm(1,:)=([q(1,4)*eye(3)-qc q(1,1:3)';-q(1,1:3) q(1,4)]*q_noise(1,:)')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid
delta=0.4;
h=1;

% Cross Product Matrix
qc=[0 -qm(1,3) qm(1,2)
   qm(1,3) 0 -qm(1,1)
   -qm(1,2) qm(1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(1,4)*eye(3)+qc;-qm(1,1:3)];
psiq=[qm(1,4)*eye(3)-qc;-qm(1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(1,4)*eye(3)+qc_d;-q_d(1,1:3)];
psiq_d=[q_d(1,4)*eye(3)-qc_d;-q_d(1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(1,:)';
werr=wm(1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(1,:)';
dq4=qm(1,:)*q_d(1,:)';
werr_hybrid(1,:)=werr';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(1,:)=torq(:)';
h_store(1)=h;

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

f1=euler_wmap_fun(xa(i,:),u(i,:),dist(i,:),intrue);
f2=euler_wmap_fun(xa(i,:)+0.5*f1'*dt,u(i,:),dist(i,:),intrue);
f3=euler_wmap_fun(xa(i,:)+0.5*f2'*dt,u(i,:),dist(i,:),intrue);
f4=euler_wmap_fun(xa(i,:)+f3'*dt,u(i,:),dist(i,:),intrue);
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion and Angular Velocity 
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Measurements
qc=[0 -q(i+1,3) q(i+1,2)
   q(i+1,3) 0 -q(i+1,1)
   -q(i+1,2) q(i+1,1) 0];
wm(i+1,:)=w(i+1,:)+w_noise(i+1,:);
qm(i+1,:)=([q(i+1,4)*eye(3)-qc q(i+1,1:3)';-q(i+1,1:3) q(i+1,4)]*q_noise(i+1,:)')';

% Cross Product Matrices 
qc=[0 -qm(i+1,3) qm(i+1,2)
   qm(i+1,3) 0 -qm(i+1,1)
   -qm(i+1,2) qm(i+1,1) 0];
qc_d=[0 -q_d(i+1,3) q_d(i+1,2)
   q_d(i+1,3) 0 -q_d(i+1,1)
   -q_d(i+1,2) q_d(i+1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(i+1,4)*eye(3)+qc;-qm(i+1,1:3)];
psiq=[qm(i+1,4)*eye(3)-qc;-qm(i+1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(i+1,4)*eye(3)+qc_d;-q_d(i+1,1:3)];
psiq_d=[q_d(i+1,4)*eye(3)-qc_d;-q_d(i+1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(i+1,:)';
werr=wm(i+1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(i+1,:)';
dq4=qm(i+1,:)*q_d(i+1,:)';
werr_hybrid(i+1,:)=werr';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(i+1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(i+1,:)=torq(:)';
torq_norm(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));
h_store(i+1)=h;

end

q_hybrid=q;
w_hybrid=w;
theta_hybrid=2*acos(abs(qm(:,1).*q_d(:,1)+qm(:,2).*q_d(:,2)+qm(:,3).*q_d(:,3)+qm(:,4).*q_d(:,4)))*180/pi;
qm_hybrid=qm;
wm_hybrid=wm;
h_store_hybrid=h_store;
u_hybrid=u;
torq_norm_hybrid=torq_norm;
h_dq4_hybrid=h_store.*qm(:,1).*q_d(:,1)+h_store.*qm(:,2).*q_d(:,2)+h_store.*qm(:,3).*q_d(:,3)+h_store.*qm(:,4).*q_d(:,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discontinuous
i1000=0;
delta=0;
h=1;

% Cross Product Matrix 
qc=[0 -qm(1,3) qm(1,2)
   qm(1,3) 0 -qm(1,1)
   -qm(1,2) qm(1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(1,4)*eye(3)+qc;-qm(1,1:3)];
psiq=[qm(1,4)*eye(3)-qc;-qm(1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(1,4)*eye(3)+qc_d;-q_d(1,1:3)];
psiq_d=[q_d(1,4)*eye(3)-qc_d;-q_d(1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(1,:)';
werr=wm(1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(1,:)';
dq4=qm(1,:)*q_d(1,:)';
werr_dis(1,:)=werr';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(1,:)=torq(:)';
h_store(1)=h;

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

f1=euler_wmap_fun(xa(i,:),u(i,:),dist(i,:),intrue);
f2=euler_wmap_fun(xa(i,:)+0.5*f1'*dt,u(i,:),dist(i,:),intrue);
f3=euler_wmap_fun(xa(i,:)+0.5*f2'*dt,u(i,:),dist(i,:),intrue);
f4=euler_wmap_fun(xa(i,:)+f3'*dt,u(i,:),dist(i,:),intrue);
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion and Angular Velocity
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Measurements
qc=[0 -q(i+1,3) q(i+1,2)
   q(i+1,3) 0 -q(i+1,1)
   -q(i+1,2) q(i+1,1) 0];
wm(i+1,:)=w(i+1,:)+w_noise(i+1,:);
qm(i+1,:)=([q(i+1,4)*eye(3)-qc q(i+1,1:3)';-q(i+1,1:3) q(i+1,4)]*q_noise(i+1,:)')';

% Cross Product Matrices 
qc=[0 -qm(i+1,3) qm(i+1,2)
   qm(i+1,3) 0 -qm(i+1,1)
   -qm(i+1,2) qm(i+1,1) 0];
qc_d=[0 -q_d(i+1,3) q_d(i+1,2)
   q_d(i+1,3) 0 -q_d(i+1,1)
   -q_d(i+1,2) q_d(i+1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(i+1,4)*eye(3)+qc;-qm(i+1,1:3)];
psiq=[qm(i+1,4)*eye(3)-qc;-qm(i+1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(i+1,4)*eye(3)+qc_d;-q_d(i+1,1:3)];
psiq_d=[q_d(i+1,4)*eye(3)-qc_d;-q_d(i+1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(i+1,:)';
werr=wm(i+1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(i+1,:)';
dq4=qm(i+1,:)*q_d(i+1,:)';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(i+1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(i+1,:)=torq(:)';
torq_norm(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));
h_store(i+1)=h;
werr_dis(i+1,:)=werr';

end

q_dis=q;
w_dis=w;
theta_dis=2*acos(abs(qm(:,1).*q_d(:,1)+qm(:,2).*q_d(:,2)+qm(:,3).*q_d(:,3)+qm(:,4).*q_d(:,4)))*180/pi;
qm_dis=qm;
wm_dis=wm;
h_store_dis=h_store;
u_dis=u;
torq_norm_dis=torq_norm;
h_dq4_dis=h_store.*qm(:,1).*q_d(:,1)+h_store.*qm(:,2).*q_d(:,2)+h_store.*qm(:,3).*q_d(:,3)+h_store.*qm(:,4).*q_d(:,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unwinding
i1000=0;
delta=3;
h=1;

% Cross Product Matrix 
qc=[0 -qm(1,3) qm(1,2)
   qm(1,3) 0 -qm(1,1)
   -qm(1,2) qm(1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(1,4)*eye(3)+qc;-qm(1,1:3)];
psiq=[qm(1,4)*eye(3)-qc;-qm(1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(1,4)*eye(3)+qc_d;-q_d(1,1:3)];
psiq_d=[q_d(1,4)*eye(3)-qc_d;-q_d(1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(1,:)';
werr=wm(1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(1,:)';
dq4=qm(1,:)*q_d(1,:)';
werr_unwind(1,:)=werr';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(1,:)=torq(:)';
h_store(1)=h;

% Main Loop
for i=1:m-1,

if (i1000==1000),
 disp(sprintf('      Program has reached point %5i',i-1))
 i1000=0;
end
i1000=i1000+1; 

f1=euler_wmap_fun(xa(i,:),u(i,:),dist(i,:),intrue);
f2=euler_wmap_fun(xa(i,:)+0.5*f1'*dt,u(i,:),dist(i,:),intrue);
f3=euler_wmap_fun(xa(i,:)+0.5*f2'*dt,u(i,:),dist(i,:),intrue);
f4=euler_wmap_fun(xa(i,:)+f3'*dt,u(i,:),dist(i,:),intrue);
xa(i+1,:)=xa(i,:)+1/6*(f1'+2*f2'+2*f3'+f4')*dt;

% Actual Quaterion and Angular Velocity
q(i+1,:)=xa(i+1,1:4);
w(i+1,:)=xa(i+1,5:7);

% Measurements
qc=[0 -q(i+1,3) q(i+1,2)
   q(i+1,3) 0 -q(i+1,1)
   -q(i+1,2) q(i+1,1) 0];
wm(i+1,:)=w(i+1,:)+w_noise(i+1,:);
qm(i+1,:)=([q(i+1,4)*eye(3)-qc q(i+1,1:3)';-q(i+1,1:3) q(i+1,4)]*q_noise(i+1,:)')';

% Cross Product Matrices 
qc=[0 -qm(i+1,3) qm(i+1,2)
   qm(i+1,3) 0 -qm(i+1,1)
   -qm(i+1,2) qm(i+1,1) 0];
qc_d=[0 -q_d(i+1,3) q_d(i+1,2)
   q_d(i+1,3) 0 -q_d(i+1,1)
   -q_d(i+1,2) q_d(i+1,1) 0];

% Xi, Psi and Attitude Matrices
xiq=[qm(i+1,4)*eye(3)+qc;-qm(i+1,1:3)];
psiq=[qm(i+1,4)*eye(3)-qc;-qm(i+1,1:3)];
attq=xiq'*psiq;
xiq_d=[q_d(i+1,4)*eye(3)+qc_d;-q_d(i+1,1:3)];
psiq_d=[q_d(i+1,4)*eye(3)-qc_d;-q_d(i+1,1:3)];
attq_d=xiq_d'*psiq_d;

% Attitude Error Matrix
att_err=attq*attq_d';

% Angular Velocity and Quaterinon Error
wd_rot=att_err*w_d(i+1,:)';
werr=wm(i+1,:)'-wd_rot;
sc=[0 -wd_rot(3) wd_rot(2)
   wd_rot(3) 0 -wd_rot(1)
   -wd_rot(2) wd_rot(1) 0];
dq=xiq_d'*qm(i+1,:)';
dq4=qm(i+1,:)*q_d(i+1,:)';

% Hysteresis
if h*dq4 >= -delta
 h=h;
else
 h=sign(dq4);
end

% Control Torque 
torq=in*att_err*w_dd(i+1,:)'+sc*in*wd_rot-kp*h*dq-kd*werr;
u(i+1,:)=torq(:)';
torq_norm(i+1)=sqrt(trapz([0:dt:t(i+1)]',u(1:i+1,1).^2+u(1:i+1,2).^2+u(1:i+1,3).^2));
h_store(i+1)=h;
werr_unwind(i+1,:)=werr';

end

q_unwind=q;
w_unwind=w;
theta_unwind=2*acos(abs(qm(:,1).*q_d(:,1)+qm(:,2).*q_d(:,2)+qm(:,3).*q_d(:,3)+qm(:,4).*q_d(:,4)))*180/pi;
qm_unwind=qm;
wm_unwind=wm;
h_store_unwind=h_store;
u_unwind=u;
torq_norm_unwind=torq_norm;
h_dq4_unwind=h_store.*qm(:,1).*q_d(:,1)+h_store.*qm(:,2).*q_d(:,2)+h_store.*qm(:,3).*q_d(:,3)+h_store.*qm(:,4).*q_d(:,4);


% Plot Results
clf
plot(t/60,h_dq4_hybrid,t/60,h_dq4_dis,'--',t/60,h_dq4_unwind,'-.')
set(gca,'fontsize',12);
axis([0 10 -1 1])
legend('Hybrid','Discontinuous','Unwinding','Location','SouthEast')
ylabel('aaa')
xlabel('Time (Min)')

pause

plot(t/60,theta_hybrid,t/60,theta_dis,'--',t/60,theta_unwind,'-.')
set(gca,'fontsize',12);
axis([0 10 0 180])
legend('Hybrid','Discontinuous','Unwinding')
ylabel('Theta')
xlabel('Time (Min)')

pause

plot(t/60,(werr_hybrid(:,1).^2+werr_hybrid(:,2).^2+werr_hybrid(:,3).^2).^(0.5),t/60,(werr_dis(:,1).^2+werr_dis(:,2).^2+werr_dis(:,3).^2).^(0.5),'--',t/60,(werr_unwind(:,1).^2+werr_unwind(:,2).^2+werr_unwind(:,3).^2).^(0.5),'-.')
set(gca,'fontsize',12);
axis([0 10 0 0.4])
legend('Hybrid','Discontinuous','Unwinding')
ylabel('dw (rad/sec)')
xlabel('Time (Min)')

pause

plot(t/60,torq_norm_hybrid,t/60,torq_norm_dis,'--',t/60,torq_norm_unwind,'-.')
set(gca,'fontsize',12);
axis([0 60 0 80])
legend('Hybrid','Discontinuous','Unwinding','Location','SouthEast')
ylabel('Torque Norm')
xlabel('Time (Min)')
