function [b_eci,b_ecef,b_ned,ti,bh,dec,dip]=mag_field(r,yr,mth,day,hr,min,sec,mln,dt,type);
% [b_eci,b_ecef,b_ned,ti,bh,dec,dip]=mag_field(r,yr,mth,day,hr,min,sec,mln,dt,type);
%
%  This program computes the reference magnetic field using the WMM model. 
%
%  The inputs are:
%       r = position vector (see type below)
%      yr = year (e.g. 2010)
%     mth = month
%     day = day
%      hr = hour (0-23)
%     min = minute (0-59)
%     sec = second (0-59)
%     mln = field order (up to 12)
%      dt = sampling interval in seconds (will be used for m > 1)
%    type = 1 for r in ECI coordinates in km
%         = 2 for r in ECEF coordinates in km
%         = 3 for r in [lat long height] in [deg deg meter]
%
%  The outputs are:
%   b_eci = magnetic field in ECI coordinates in nT (m x 3)
%  b_ecef = magnetic field in ECEF coordinates in nT (m x 3)
%   b_ned = magnetic field in NED coordinates in nT (m x 3)
%      ti = total field intensity in nT (m x 1)
%      bh = horizontal field intensity in nT (m x 1)
%     dec = declination in deg (m x 1)
%     dip = inclination in deg (m x 1)

% John L. Crassidis 5/28/06

% Note: this program was converted from FORTRAN code provided by the
% National Geophysical Data Center. See http://www.ngdc.noaa.gov/seg/WMM/.

% Epoch
epoch=2010.00;

% Load Data for 2010
wmm=wmm2010_data;

% Constants for Degrees/Radians Conversions
dtr=pi/180;
rtd=180/pi;

% Get Length of Data
[mout,mout3]=size(r);
if mout3~=3, disp('r must be (m x 3)');end

% Initialize Field Vectors
bt_store=zeros(mout,1);br_store=zeros(mout,1);bp_store=zeros(mout,1);
t=[0:dt:(mout-1)*dt]';

% Convert from ECI to Lat, Long and Height if Needed
if type == 1
 r_ecef=eci2ecef(r,yr,mth,day,hr,min,sec+t);
 [lat,long,height]=ecef2llh(r_ecef);height=height/1000;
end

% Convert from ECEF to Lat, Long and Height if Needed
if type == 2
 [lat,long,height]=ecef2llh(r);height=height/1000;
end

% Get Lat, Long and Height
if type == 3
 lat=r(:,1);
 long=r(:,2);
 height=r(:,3)/1000;
end

% Convert to Radians
rlon=long*dtr;
rlat=lat*dtr;

% Get Sidereal Time and Ascending Node
theta=sidereal(yr,mth,day,hr,min,sec+t)*dtr;
alpha=rlon+theta;
c_alpha=cos(alpha);s_alpha=sin(alpha);

% Leap Year
if (mod(yr,400)&&~mod(yr,100))
% leapyear = false;
 ndays = 365;
elseif ~mod(yr,4)
% leapyear = true;
 ndays = 366;
else
% leapyear = false;
 ndays = 365;
end

% Get Days Past First of Year
o1=ones(mout,1);
day_of_year = datenum(yr*o1,mth*o1,day*o1,hr*o1,min*o1,sec+t)-datenum(yr,1,1);

% Find Days That Go Into Next Year
i_days=find(day_of_year > ndays);

% Get Decimal Year Time
if isempty(i_days) == 1
 time = yr + day_of_year/ndays;
else
 % Check Leap Year of New Year
 time = yr + day_of_year/ndays;
 yr_new=floor(time(i_days(1)));
 if (mod(yr_new,400)&&~mod(yr_new,100))
  % leapyear = false; 
  ndays_new = 365;
 elseif ~mod(yr_new,4)
  % leapyear = true;
  ndays_new = 366;
 else
  % leapyear = false;
  ndays_new = 365;
 end 
 time=zeros(mout,1);
 time(1:i_days(1)-1) = yr + day_of_year(1:i_days(1)-1)/ndays;
 time(i_days(1):mout) = yr + (day_of_year(i_days:mout)-ndays+ndays_new)/ndays_new;
end

% Get Spherical Harmonic Coefficients from Data
c=zeros(13,13);cd=zeros(13,13);
for i=1:length(wmm)
 n=wmm(i,1);
 m=wmm(i,2);
 if (m <= n)
  c(n+1,m+1)=wmm(i,3);
  cd(n+1,m+1)=wmm(i,5);
 end
 if (m ~= 0)
  c(m,n+1)=wmm(i,4);
  cd(m,n+1)=wmm(i,6);
 end
end

% Convert Schmidt Normalized Gauss Coefficients to Unnormalized
snorm=zeros(mln+1,mln+1);
k=zeros(mln+1,mln+1);
fn=zeros(mln+1,1);
fm=zeros(mln+1,1);

snorm(1,1)=1;
for n=2:mln+1
 snorm(n,1)=snorm(n-1,1)*(2*(n-1)-1)/(n-1);
 j=2;
 for m=1:n
  k(n,m)=((n-2)^2-(m-1)^2)/((2*(n-1)-1)*(2*(n-1)-3));
  if m > 1
   flnmj=((n-m+1)*j)/(n+m-2);
   snorm(n,m)=snorm(n,m-1)*sqrt(flnmj);
   j=1;
   c(m-1,n)=snorm(n,m)*c(m-1,n);
   cd(m-1,n)=snorm(n,m)*cd(m-1,n);
  end
  c(n,m)=snorm(n,m)*c(n,m);
  cd(n,m)=snorm(n,m)*cd(n,m);  
 end
 fn(n)=n;
 fm(n)=n-1;
end
k(2,2)=0;

% Constant for WGS-84
a=6378.137;
b=6356.7523142;
re=6371.2;
a2=a^2;
b2=b^2;
c2=a2-b2;
a4=a2^2;
b4=b2^2;
c4=a4-b4;

% Sine and Cosine of Lat and Long
srlon=sin(rlon);
srlat=sin(rlat);
crlon=cos(rlon);
crlat=cos(rlat);
srlon2=srlon.^2;
srlat2=srlat.^2;
crlon2=crlon.^2;
crlat2=crlat.^2;

% Convert from Geodetic to Spherical Coordinates
q=(a2-c2*srlat2).^(0.5);
q1=height.*q;
q2=((q1+a2)./(q1+b2)).^2;
ct=srlat./((q2.*crlat2+srlat2).^(0.5));
st=(1-ct.^2).^(0.5);
r2=height.^2+2*q1+(a4-c4*srlat2)./(q.^2);
r=(r2).^(0.5);
d=(a2*crlat2+b2*srlat2).^(0.5);
ca=(height+d)./r;
sa=c2*crlat.*srlat./(r.*d);

% Main Loop
for i = 1:mout

% Update Time
dt_change=time(i)-epoch;

% Initialize Quantities
sp=zeros(13,1);
cp=zeros(13,1);
p=zeros(13,13);
dp=zeros(13,13);
pp=zeros(13,1);

p(1,1)=1;
pp(1)=1;
dp(1,1)=0;

sp(1)=0;
cp(1)=1;
sp(2)=srlon(i);
cp(2)=crlon(i);

for m=3:mln+1;
 sp(m)=sp(2)*cp(m-1)+cp(2)*sp(m-1);
 cp(m)=cp(2)*cp(m-1)-sp(2)*sp(m-1);
end

aor=re/r(i);
ar=aor^2;

br=0;
bt=0;
bp=0;
bpp=0;

tc=zeros(mln+1,mln+1);

% Loop for n and m
for n=2:mln+1
 ar=ar*aor;
 for m=1:n
     
% Compute Unnormalized Associated Legendre Polynomials and Derivatives
  if (n == m)
   p(n,m)=st(i)*p(n-1,m-1);
   dp(n,m)=st(i)*dp(n-1,m-1)+ct(i)*p(n-1,m-1);
  elseif (n == 2) & (m == 1)
   p(n,m)=ct(i)*p(n-1,m);
   dp(n,m)=ct(i)*dp(n-1,m)-st(i)*p(n-1,m);
  elseif (n > 2) & (n ~= m)
   if (m > n-2) 
    p(n-2,m)=0;
    dp(n-2,m)=0;
   end
   p(n,m)=ct(i)*p(n-1,m)-k(n,m)*p(n-2,m);
   dp(n,m)=ct(i)*dp(n-1,m)-st(i)*p(n-1,m)-k(n,m)*dp(n-2,m);
  end
  
% Time Adjust the Gauss Coefficients
  tc(n,m)=c(n,m)+dt_change*cd(n,m);   
  if (m ~= 1)
   tc(m-1,n)=c(m-1,n)+dt_change*cd(m-1,n);
  end
  
% Accumulate Terms of the Spherical Harmonic Expansions
  par=ar*p(n,m);
  if (m == 1)
   temp1=tc(n,m)*cp(m);
   temp2=tc(n,m)*sp(m);
  else
   temp1=tc(n,m)*cp(m)+tc(m-1,n)*sp(m);
   temp2=tc(n,m)*sp(m)-tc(m-1,n)*cp(m);
  end
  bt=bt-ar*temp1*dp(n,m);
  bp=bp+fm(m)*temp2*par;
  br=br+fn(n)*temp1*par; 
 
% Special Case: North/South Geographic Poles
  if (st(i) == 0) & (m == 2)
   if (n == 2)
    pp(n)=pp(n-1);
   else
    pp(n)=ct(i)*pp(n-1)-k(n,m)*pp(n-2);
   end
   parp=ar*pp(n);
   bpp=bpp+fm(m)*temp2*parp;
  end
 end
end
 
if (st(i) == 0)
 bp=bpp;
else
 bp=bp/st(i);
end

% Store Field Variables
bt_store(i)=bt;
br_store(i)=br;
bp_store(i)=bp;

end

% Rotate from Spherical to Geodetic Coordinates
b_ned=[-bt_store.*ca-br_store.*sa bp_store bt_store.*sa-br_store.*ca];

% Get ECI Coordinates
brbt=(br_store.*st+bt_store.*ct);
b_eci=[brbt.*c_alpha-bp_store.*s_alpha brbt.*s_alpha+bp_store.*c_alpha br_store.*ct-bt_store.*st];

% Horizontal, Total Field, Declination and Dip
bh=(b_ned(:,1).^2+b_ned(:,2).^2).^(0.5);
ti=(bh.^2+b_ned(:,3).^2).^(0.5);
dec=atan2(b_ned(:,2),b_ned(:,1))*rtd;
dip=atan2(b_ned(:,3),bh)*rtd;

% Get ECEF Coordinates
b_ecef=eci2ecef(b_eci,yr,mth,day,hr,min,sec+t);

return

function out=wmm2010_data

out=[1  0  -29496.6       0.0       11.6        0.0
  1  1   -1586.3    4944.4       16.5      -25.9
  2  0   -2396.6       0.0      -12.1        0.0
  2  1    3026.1   -2707.7       -4.4      -22.5
  2  2    1668.6    -576.1        1.9      -11.8
  3  0    1340.1       0.0        0.4        0.0
  3  1   -2326.2    -160.2       -4.1        7.3
  3  2    1231.9     251.9       -2.9       -3.9
  3  3     634.0    -536.6       -7.7       -2.6
  4  0     912.6       0.0       -1.8        0.0
  4  1     808.9     286.4        2.3        1.1
  4  2     166.7    -211.2       -8.7        2.7
  4  3    -357.1     164.3        4.6        3.9
  4  4      89.4    -309.1       -2.1       -0.8
  5  0    -230.9       0.0       -1.0        0.0
  5  1     357.2      44.6        0.6        0.4
  5  2     200.3     188.9       -1.8        1.8
  5  3    -141.1    -118.2       -1.0        1.2
  5  4    -163.0       0.0        0.9        4.0
  5  5      -7.8     100.9        1.0       -0.6
  6  0      72.8       0.0       -0.2        0.0
  6  1      68.6     -20.8       -0.2       -0.2
  6  2      76.0      44.1       -0.1       -2.1
  6  3    -141.4      61.5        2.0       -0.4
  6  4     -22.8     -66.3       -1.7       -0.6
  6  5      13.2       3.1       -0.3        0.5
  6  6     -77.9      55.0        1.7        0.9
  7  0      80.5       0.0        0.1        0.0
  7  1     -75.1     -57.9       -0.1        0.7
  7  2      -4.7     -21.1       -0.6        0.3
  7  3      45.3       6.5        1.3       -0.1
  7  4      13.9      24.9        0.4       -0.1
  7  5      10.4       7.0        0.3       -0.8
  7  6       1.7     -27.7       -0.7       -0.3
  7  7       4.9      -3.3        0.6        0.3
  8  0      24.4       0.0       -0.1        0.0
  8  1       8.1      11.0        0.1       -0.1
  8  2     -14.5     -20.0       -0.6        0.2
  8  3      -5.6      11.9        0.2        0.4
  8  4     -19.3     -17.4       -0.2        0.4
  8  5      11.5      16.7        0.3        0.1
  8  6      10.9       7.0        0.3       -0.1
  8  7     -14.1     -10.8       -0.6        0.4
  8  8      -3.7       1.7        0.2        0.3
  9  0       5.4       0.0        0.0        0.0
  9  1       9.4     -20.5       -0.1        0.0
  9  2       3.4      11.5        0.0       -0.2
  9  3      -5.2      12.8        0.3        0.0
  9  4       3.1      -7.2       -0.4       -0.1
  9  5     -12.4      -7.4       -0.3        0.1
  9  6      -0.7       8.0        0.1        0.0
  9  7       8.4       2.1       -0.1       -0.2
  9  8      -8.5      -6.1       -0.4        0.3
  9  9     -10.1       7.0       -0.2        0.2
 10  0      -2.0       0.0        0.0        0.0
 10  1      -6.3       2.8        0.0        0.1
 10  2       0.9      -0.1       -0.1       -0.1
 10  3      -1.1       4.7        0.2        0.0
 10  4      -0.2       4.4        0.0       -0.1
 10  5       2.5      -7.2       -0.1       -0.1
 10  6      -0.3      -1.0       -0.2        0.0
 10  7       2.2      -3.9        0.0       -0.1
 10  8       3.1      -2.0       -0.1       -0.2
 10  9      -1.0      -2.0       -0.2        0.0
 10 10      -2.8      -8.3       -0.2       -0.1
 11  0       3.0       0.0        0.0        0.0
 11  1      -1.5       0.2        0.0        0.0
 11  2      -2.1       1.7        0.0        0.1
 11  3       1.7      -0.6        0.1        0.0
 11  4      -0.5      -1.8        0.0        0.1
 11  5       0.5       0.9        0.0        0.0
 11  6      -0.8      -0.4        0.0        0.1
 11  7       0.4      -2.5        0.0        0.0
 11  8       1.8      -1.3        0.0       -0.1
 11  9       0.1      -2.1        0.0       -0.1
 11 10       0.7      -1.9       -0.1        0.0
 11 11       3.8      -1.8        0.0       -0.1
 12  0      -2.2       0.0        0.0        0.0
 12  1      -0.2      -0.9        0.0        0.0
 12  2       0.3       0.3        0.1        0.0
 12  3       1.0       2.1        0.1        0.0
 12  4      -0.6      -2.5       -0.1        0.0
 12  5       0.9       0.5        0.0        0.0
 12  6      -0.1       0.6        0.0        0.1
 12  7       0.5       0.0        0.0        0.0
 12  8      -0.4       0.1        0.0        0.0
 12  9      -0.4       0.3        0.0        0.0
 12 10       0.2      -0.9        0.0        0.0
 12 11      -0.8      -0.2       -0.1        0.0
 12 12       0.0       0.9        0.1        0.0];

function [lat,long,height]=ecef2llh(r_ecef)
% [lat,long,height] = ecef2llh(r_ecef)
%
%  This program converts ECEF quantities into long, lat and height.
%
%  The input is:
%   r_ecef = ECEF quantities in km (m x 3)
%
%  The outputs are:
%      lat = geodetic latitude in degrees (m x 1)
%     long = longitude in degrees (m x 1)
%   height = height in meters (m x 1)

% John L. Crassidis 6/7/06

% This program uses a closed-form solution, which is valid as long as
% does not fall within 43 km of the Earth center. Details can be found at:
% Zhu, J., "Exact Conversion of Earth-Centered-Earth-Fixed Coordinates to 
% Geodetic Coordinates," Journal of Guidance, Control and Dynamics, 
% Vol. 16, No. 2, March-April 1993, pp. 389-391.

% Check Size of Input Vector
rad2deg=180/pi;
[mout,mout1]=size(r_ecef);
if mout1 ~= 3
 disp(' r_ecef must be m x 3, where m is the total number of points')
end

% Constant for WGS-84
a=6378.137;
b=6356.7523142;
e2=0.00669437999;

% Conversion Parameters
w=(r_ecef(:,1).^2+r_ecef(:,2).^2).^(0.5);
l=e2/2;m=(w/a).^2;n=((1-e2)*r_ecef(:,3)/b).^2;
i=-(2*l^2+m+n)/2;k=l^2*(l^2-m-n);
mnl2=m.*n*l^2;
q=(m+n-4*l^2).^3/216+mnl2;
d=((2*q-mnl2).*mnl2).^(0.5);
beta=i/3-(q+d).^(1/3)-(q-d).^(1/3);
t=((beta.^2-k).^(0.5)-(beta+i)/2).^(0.5)-sign(m-n).*((beta-i)/2).^(0.5);
w1=w./(t+l);z1=(1-e2)*r_ecef(:,3)./(t-l);
 
% Compute Longitude
jj0=find((w-r_ecef(:,1))==0);
jj=find((w-r_ecef(:,1))~=0);
long=zeros(mout,1);
if isempty(jj) == 0
 long(jj,:)=2*atan((w(jj)-r_ecef(jj,1))./r_ecef(jj,2))*rad2deg;
end
if isempty(jj0) == 0
 long(jj0,:)=zeros(length(jj0),1);
end
 
% Compute Latitude
j0=find(w==0);
j=find(w~=0);
lat=zeros(mout,1);height=zeros(mout,1);
if isempty(j) == 0
 lat(j,:)=atan(z1(j)./((1-e2)*w1(j)))*rad2deg;
 height(j,:)=sign(t(j)-1+l).*((w(j)-w1(j)).^2+(r_ecef(j,3)-z1(j)).^2).^(0.5)*1000;
end
if isempty(j0) == 0
 lat(j0)=sign(r_ecef(j0,3))*90;
 h(j0)=sign(r_ecef(j0,3))*r_ecef(j0,3)-b;
end

return

function r_ecef=llh2ecef(lat,long,height)
% r_ecef = llh2ecef(lat,long,height)
%
%  This program converts long, lat and height into ECEF quantities.
%
%  The inputs are:
%      lat = geodetic latitude in degrees (m x 1)
%     long = longitude in degrees (m x 1)
%   height = height in meters (m x 1)
%
%  The output is:
%   r_ecef = ECEF quantities in km (m x 3)

% John L. Crassidis 6/7/06

deg2rad=pi/180;

% Get Length of Data
mout=length(lat);
if length(long) ~= mout
 disp(' Lat, Long and Height Must Have Same Length')
end
if length(height) ~= mout
 disp(' Lat, Long and Height Must Have Same Length')
end

% Constant for WGS-84
a=6378.137;
b=6356.7523142;
e2=0.00669437999;

% Convert Lat and Long to Degrees and Height to Km
lat=lat(:)*deg2rad;
long=long(:)*deg2rad;
height=height(:)/1000;

% Compute ECEF Quantities
n=a./((1-e2*sin(lat).^2).^(0.5));
x=(n+height).*cos(lat).*cos(long);
y=(n+height).*cos(lat).*sin(long);
z=(n*(1-e2)+height).*sin(lat);
r_ecef=[x y z];

return

function r_eci=ecef2eci(r_ecef,yr,mth,day,hr,min,sec)
% r_eci=ecef2eci(r_ecef,yr,mth,day,hr,min,sec)
%
%  This program converts ECEF into ECI quantities.
%
%  The inputs are:
%   r_ecef = ECEF quantities (m x 3)
%       yr = year (e.g. 2005), can be (m x 1)
%      mth = month, can be (m x 1)
%      day = day, can be (m x 1)
%       hr = hour (0-23), can be (m x 1)
%      min = minute (0-59),can be (m x 1)
%      sec = second (0-59), can be (m x 1)
%
%  The output is:
%     r_eci= ECEF quantities (m x 3)

% John L. Crassidis 6/7/06

% Sidereal Time is Computed using Eq. (12.4) from Meeus, J., 
% Astronomical Algorithms, 2nd Edition, Willmann-Bell, Inc., 
% Richmond, VA, 1998.

% Check Size of Input Vector
[mout,mout1]=size(r_ecef);
if mout1 ~= 3
 disp(' r_ecef must be m x 3, where m is the total number of points')
end

% Compute Day-of-Year
doy = day + 31*(mth-1) - fix(2.2+0.4*mth).*(mth>2) + (mth>2).*(~rem(yr,4));

% Compute Julian Date
jd = 2415020 + (yr-1900)*365 + fix((yr-1901)/4) + (doy-1) + 0.5 + ...
       (3600*hr + 60*min + sec) / 86400;

% Compute T
t=(jd-2451545)/36525;

% Compute Sidereal Angle
theta=(280.46061837+360.98564736629*(jd-2451545)+0.000387933*t.^2-t.^3/38710000)*pi/180;

% Compute ECI Components
cost=cos(theta);sint=sin(theta);
r_eci(:,1)=cost.*r_ecef(:,1)-sint.*r_ecef(:,2);
r_eci(:,2)=sint.*r_ecef(:,1)+cost.*r_ecef(:,2);
r_eci(:,3)=r_ecef(:,3);

return

function r_ecef=eci2ecef(r_eci,yr,mth,day,hr,min,sec)
% r_ecef=eci2ecef(r_eci,yr,mth,day,hr,min,sec)
%
%  This program converts ECI into ECEF quantities.
%
%  The inputs are:
%     r_eci= ECEF quantities (m x 3)
%       yr = year (e.g. 2005), can be (m x 1)
%      mth = month, can be (m x 1)
%      day = day, can be (m x 1)
%       hr = hour (0-23), can be (m x 1)
%      min = minute (0-59),can be (m x 1)
%      sec = second (0-59), can be (m x 1)
%
%  The output is:
%   r_ecef = ECEF quantities (m x 3)

% John L. Crassidis 6/7/06

% Sidereal Time is Computed using Eq. (12.4) from Meeus, J., 
% Astronomical Algorithms, 2nd Edition, Willmann-Bell, Inc., 
% Richmond, VA, 1998.

% Check Size of Input Vector
[mout,mout1]=size(r_eci);
if mout1 ~= 3
 disp(' r_eci must be m x 3, where m is the total number of points')
end

% Compute Day-of-Year
doy = day + 31*(mth-1) - fix(2.2+0.4*mth).*(mth>2) + (mth>2).*(~rem(yr,4));

% Compute Julian Date
jd = 2415020 + (yr-1900)*365 + fix((yr-1901)/4) + (doy-1) + 0.5 + ...
       (3600*hr + 60*min + sec) / 86400;

% Compute T
t=(jd-2451545)/36525;

% Compute Sidereal Angle
theta=(280.46061837+360.98564736629*(jd-2451545)+0.000387933*t.^2-t.^3/38710000)*pi/180;

% Compute ECEF Components
cost=cos(theta);sint=sin(theta);
r_ecef(:,1)=cost.*r_eci(:,1)+sint.*r_eci(:,2);
r_ecef(:,2)=-sint.*r_eci(:,1)+cost.*r_eci(:,2);
r_ecef(:,3)=r_eci(:,3);

return

function theta=sidereal(yr,mth,day,hr,min,sec)
% theta = sidereal(yr,mth,day,hr,min,sec)
%
%  This program computes the mean sidereal time in degrees.
%
%  The inputs are:
%       yr = year (e.g. 2005), can be (m x 1)
%      mth = month, can be (m x 1)
%      day = day, can be (m x 1)
%       hr = hour (0-23), can be (m x 1)
%      min = minute (0-59),can be (m x 1)
%      sec = second (0-59), can be (m x 1)
%
%  The output is:
%    theta = sidereal time in degrees (m x 3)

% John L. Crassidis 6/7/06

% Sidereal Time is Computed using Eq. (12.4) from Meeus, J., 
% Astronomical Algorithms, 2nd Edition, Willmann-Bell, Inc., 
% Richmond, VA, 1998.


% Compute Day-of-Year
doy = day + 31*(mth-1) - fix(2.2+0.4*mth).*(mth>2) + (mth>2).*(~rem(yr,4));

% Compute Julian Date
jd = 2415020 + (yr-1900)*365 + fix((yr-1901)/4) + (doy-1) + 0.5 + ...
       (3600*hr + 60*min + sec) / 86400;

% Compute T
t=(jd-2451545)/36525;

% Compute Sidereal Angle
theta=280.46061837+360.98564736629*(jd-2451545)+0.000387933*t.^2-t.^3/38710000;
mout=length(theta);

% Make Sidereal Angle from -360 to 360 Degrees
if (theta(1)<-360)
 theta=theta+abs(fix(theta(1)/360))*360;
end
if (theta(1)>360)
 theta=theta-abs(fix(theta(1)/360))*360;
end

return