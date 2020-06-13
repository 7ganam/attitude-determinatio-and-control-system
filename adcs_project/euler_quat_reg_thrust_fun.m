function f = euler_quat_reg_thrust_fun(x,j,invj,u);
q=x(1:4);q=q(:);w=x(5:7);w=w(:);

% Function
wc=[0 -w(3) w(2)
   w(3) 0 -w(1)
  -w(2) w(1) 0];
om=[-wc w;-w' 0];
f=[0.5*om*q;-invj*wc*j*w+invj*u(:)];