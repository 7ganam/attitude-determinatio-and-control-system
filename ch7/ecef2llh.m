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
 long(jj,:)=2*atan2((w(jj)-r_ecef(jj,1)),r_ecef(jj,2))*rad2deg;
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