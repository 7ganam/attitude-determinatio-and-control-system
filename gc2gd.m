function [long,lat,alt,gst,r_ecef]=gc2gd(r,yr,mth,day,hr,min,sec,dt,tf);
%function [long,lat,alt,alp,gst,r_ecef]=gc2gd(r,yr,mth,day,hr,min,sec,dt,tf);
%
% This function converts from Geocentric to Geodetic quantities.
%
%  The inputs are:
%       r = Geocentric coodinates (km) (m x 3)
%      yr = year, e.g. 1995
%     mth = month, e.g. Jan=1, Feb=2, etc.
%     day = day, e.g. 1-31
%      hr = hour, e.g. 0-23
%     min = minutes, e.g. 0-59
%     sec = seconds, e.g. 0-59
%      dt = sampling interval (sec)
%      tf = run time (sec)
%
%  The outputs are:
%    long = longitude (deg) (between -180 W and 180 E) (m x 1)
%     lat = latitude (deg) (m x 1)
%     alt = Geodetic altitude (km) (m x 1)
%     gst = Sidereal time (deg) (m x 1)
%  r_ecef = ECEF coordinates (km) (m x 3)

% John L. Crassidis 5/17/11

% Time vector
t=[0:dt:tf]';

% Sideral time at Greenwich from Vallado (third edition), p. 191
jdate=julian(yr,mth,day,hr,min,sec+t);
tdays=jdate-2451545;
jcent=tdays/36525;
%ut=((sec+t)/60/60+min/60+hr)*360/24;
%gst=99.6910+36000.7689*jcent+0.0004*jcent.*jcent+ut;
% Vallado p. 194-195
gst_sec=67310.54841+(876600*3600+8640184.812866)*jcent+0.093104*jcent.^2-6.2e-6*jcent.^3;
gst=rem(gst_sec,86400)/240;
i_neg=find(gst<0);gst(i_neg)=gst(i_neg)+360;

theta=gst*pi/180;
r_ecef=[r(:,1).*cos(theta)+r(:,2).*sin(theta) -r(:,1).*sin(theta)+r(:,2).*cos(theta) r(:,3)];

% Sofair, I., "Improved Method for Calculating Exact Geodetic Latitude
% and Altitude Revisited," Journal of Guidance, Control, and Dynamics,
% Vol. 23, No. 2, 2000, p. 369-369

a=6378137/1000;
b=6356752.3142/1000;
e2=1-b^2/a^2;
eps2=a^2/b^2-1;

r0=(r_ecef(:,1).^2+r_ecef(:,2).^2).^(0.5);
p=abs(r_ecef(:,3))/eps2;
s=r0.^2/(e2*eps2);
q=p.^2-b^2+s;

u=p./(q.^(0.5));
v=b^2*u.^2./q;
bigp=27*v.*s./q;

bigq=((bigp+1).^(0.5)+bigp.^(0.5)).^(2/3);
t=(1+bigq+1./bigq)/6;
c=(u.^2-1+2*t).^(0.5);
w=(c-u)/2;
z=sign(r_ecef(:,3)).*(q.^(0.5)).*(w+((t.^2+v).^(0.5)-u.*w-t/2-1/4).^(0.5));
ne=a.*(1+eps2*z.^2/b^2).^(0.5);

lat=asin((eps2+1)*(z./ne));
alt=r0.*cos(lat)+r_ecef(:,3).*sin(lat)-(1./ne)*a^2;
lat=lat*180/pi;

long=atan2(r_ecef(:,2),r_ecef(:,1))*180/pi;