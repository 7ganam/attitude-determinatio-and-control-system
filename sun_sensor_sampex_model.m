function sun_bm=sun_sensor_sampex_model(sun_b,n,h,res,sig_noise)

d2r=pi/180;
sig_noise_rad=sig_noise*d2r;
res_rad=res*d2r;

% z and x Displacement Measurements (equal to a and b in Wertz, pg 224)
z=h*sun_b(3)/sqrt(n^2-sun_b(1)^2-sun_b(3)^2);
s3=-n*z/sqrt((z*sun_b(1)/sun_b(3))^2+z^2+h^2);
z=-z*sign(s3*sun_b(3));
x=z*sun_b(1)/sun_b(3);

% True Azimuth and Coelevation
phi=atan2(x,z);
theta=atan(n*sqrt(x^2+z^2)/sqrt(h^2-(n^2-1)*(x^2+z^2)));

% Measurements
phim0=phi+sig_noise_rad*randn(1);
phim=round(phim0/res_rad)*res_rad;
thetam0=theta+sig_noise_rad*randn(1);
thetam=round(thetam0/res_rad)*res_rad;
tan_thetam=tan(thetam);
tan_betam=tan_thetam*cos(phim);
tan_alpham=tan_thetam*sin(phim);

% Get Unit Vector
sun_bm=1/sqrt(1+tan_alpham^2+tan_betam^2)*[tan_alpham;1;tan_betam];

% Get Correct Sign
sun_bm(2)=sun_bm(2)*sign(sun_bm(2)*sun_b(2));


