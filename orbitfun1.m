function f=orbitfun1(x,mu,rearth);

% Written by John L. Crassidis 5/11

% Function Routine for Orbital Equations (includes J2)
% Note this program assumes position in km and velocity in km/sec.
[m,n]=size(x);
f=zeros(m,n);

% Compute Magnitudes
r2=(x(1,:).^2+x(2,:).^2+x(3,:).^2);
r=r2.^(1/2);
r_3_2=r.^3;

% J2 Effect and Acceleration Coefficient
j2=1082.63e-6;
z_r=x(3,:)./r;
z_r2=z_r.^2;
j2_coeff=3/2*j2*(mu./r2).*(rearth^2./r2);

% Function
f(1,:)=x(4,:);
f(2,:)=x(5,:);
f(3,:)=x(6,:);
f(4,:)=-mu./r_3_2.*x(1,:)-j2_coeff.*(1-5*z_r2).*x(1,:)./r;
f(5,:)=-mu./r_3_2.*x(2,:)-j2_coeff.*(1-5*z_r2).*x(2,:)./r;
f(6,:)=-mu./r_3_2.*x(3,:)-j2_coeff.*(3-5*z_r2).*z_r;