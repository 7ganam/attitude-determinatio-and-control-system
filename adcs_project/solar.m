function [sun_i,ra_dec]=solar(yr,mth,day,hr,min,sec,dt,tf);
%function [sun_i,ra_dec]=solar(yr,mth,day,hr,min,sec,dt,tf);
%
% This function computes the solar ephemeris given date and time
%
%  The inputs are:      
%         yr = year, e.g. 1995
%        mth = month, e.g. Jan=1, Feb=2, etc.
%        day = day, e.g. 1-31
%         hr = hour, e.g. 0-23
%        min = minutes, e.g. 0-59
%        sec = seconds, e.g. 0-59
%         dt = sampling interval (sec)
%         tf = run time (sec)
%  The outputs are:
%      sun_i = x,y,z components of ECI unit vector
%     ra_dec = right ascension and declination (deg)
%

% John L. Crassidis 4/24/95

t=[0:dt:tf]';
sec=sec+t;
d2r = pi/180;
r2d = 1/d2r;

% Days past 1990
days = [ 1 32 60 91 121 152 182 213 244 274 305 335 ]';
day90 = days(mth) + day - 1 + 365*(yr - 1990);
if (day90 >= 790)
  day90 = day90 + 1;
end
day90 = day90 + (3600*hr + 60*min + sec)/86400;

% Number of days from J2000
n = -3653.5 + day90;		

% Mean longitude (deg)		
l = 280.460 + 0.9856474*n;			
while l < 0
  l = l + 360;
end

% Mean anomaly (deg)
g = 357.528 + 0.9856003*n;			
while g < 0
  g = g + 360;
end

% Ecliptic longitude (deg)
long = l + 1.915*sin(d2r*g) + 0.020*sin(2*d2r*g);	

% Obliquity of ecliptic (deg)
e = 23.439 - 0.0000004*n;			

% Right ascension (deg)
ra = r2d*atan(cos(d2r*e).*tan(d2r*long));		
while ra < 0
  ra = ra + 360;
end

% Declination (deg)
dec = r2d*asin(sin(d2r*e).*sin(d2r*long));		

% X,Y,Z components of ECI vector
sx = cos(d2r*long);				
sy = cos(d2r*e).*sin(d2r*long);			
sz = sin(d2r*e).*sin(d2r*long);

% Output
sun_i=[sx sy sz];
ra_dec=[ra dec];
