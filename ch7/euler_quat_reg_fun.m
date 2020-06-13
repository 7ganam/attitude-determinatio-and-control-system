function f = euler_quat_reg_fun(x,j,invj,qd,kp,kd);
q=x(1:4);
q=q(:);

% Error Quaternion
quat_mult_mat=[qd(4) qd(3) -qd(2) -qd(1)
              -qd(3) qd(4) qd(1) -qd(2)
               qd(2) -qd(1) qd(4) -qd(3)
               qd(1) qd(2) qd(3) qd(4)];
qerr=quat_mult_mat*q;

% Torque Command
w=x(5:7);
w=w(:);
%u=-kp*qerr(1:3)-kd*w;
u=-kp*sign(qerr(4))*qerr(1:3)-kd*w;
% Try plus or minus
%u=-kp*sign(qerr(4))*qerr(1:3)-kd*(1+qerr(1:3)'*qerr(1:3))*w;
%u=-kp*sign(qerr(4))*qerr(1:3)-kd*(1-qerr(1:3)'*qerr(1:3))*w;

% Function
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
  -w(2) w(1) 0];

om=[-wc w;-w' 0];

f=[0.5*om*q; -invj*wc*j*w+invj*u(:)];