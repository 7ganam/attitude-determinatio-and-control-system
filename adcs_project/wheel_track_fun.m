function [f,torq,ss,q_d,w_d]=wheel_track_fun(x,t,theta_d,phir_d,psir_d,d,j,jtrue,k,g,eps)

% Pre-Allocate Space
f=zeros(10,1);d=d(:);q_d=zeros(4,1);w_d=zeros(3,1);w_dd=zeros(3,1);

% Quaternion, Angular Velocity and Wheel Speed
q=[x(1);x(2);x(3);x(4)];
w=[x(5);x(6);x(7)];
wh=[x(8);x(9);x(10)];

% Desired Roll and Yaw
phi_d=phir_d*t;
psi_d=psir_d*t;

% Desired Quaternion 
q_d(1)=sin(theta_d/2)*cos((phi_d-psi_d)/2);
q_d(2)=sin(theta_d/2)*sin((phi_d-psi_d)/2);
q_d(3)=cos(theta_d/2)*sin((phi_d+psi_d)/2);
q_d(4)=cos(theta_d/2)*cos((phi_d+psi_d)/2);

% Desired Angular Velocity and Derivative
w_d(1)=sin(theta_d)*sin(psi_d)*phir_d;
w_d(2)=sin(theta_d)*cos(psi_d)*phir_d;
w_d(3)=cos(theta_d)*phir_d+psir_d;

w_dd(1)=psir_d*sin(theta_d)*cos(psi_d)*phir_d;
w_dd(2)=-psir_d*sin(theta_d)*sin(psi_d)*phir_d;
w_dd(3)=0;

% Cross Product Matrices 
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
   -w(2) w(1) 0];
qc=[0 -q(3) q(2)
   q(3) 0 -q(1)
   -q(2) q(1) 0];
qc_d=[0 -q_d(3) q_d(2)
   q_d(3) 0 -q_d(1)
   -q_d(2) q_d(1) 0];

% Xi Matrices
xiq=[q(4)*eye(3)+qc;-q(1:3)'];
xiq_d=[q_d(4)*eye(3)+qc_d;-q_d(1:3)'];

% Sliding Surface
ss=w-w_d+k*sign(q_d'*q)*xiq_d'*q;

% Control Torque for Wheel
usat=[sat1(ss(1),eps);sat1(ss(2),eps);sat1(ss(3),eps)];
dqd=xiq_d'*xiq*w-xiq'*xiq_d*w_d;
torq=wc*(j*w+wh)-j*(0.5*k*sign(q_d'*q)*dqd-w_dd+g*usat);


% Integration Functions
om=[-wc w;-w' 0];
f(1:4,:)=0.5*om*q;
f(5:7,:)=inv(jtrue)*(-wc*jtrue*w+torq);
f(8:10,:)=-wc*wh-torq;