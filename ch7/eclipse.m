function out=eclipse(r,sun_b);
%function out=eclipse(r,sun_b);
%
% This function determines eclipse conditions.
%
%  The inputs are:
%         r = spacecraft position (km)
%     sun_b = inertial sun vector (assumed normalized)
%
%  The output is:
%       out = 1 for sun available
%           = 0 for eclipse

% John L. Crassidis 4/24/95

% Earth and Sun radii constants
rs=6.9599e5;
rp=6378.140;

% Find magnitude and normalize
dsmag=149.6e6;
dpmag=vecnorm(r);

% Eclipse parameters
rhos=asin(rs./dsmag);
rhop=asin(rp./dpmag);
theta=acos(sun_b(:,1).*dpmag+sun_b(:,2).*dpmag+sun_b(:,3).*dpmag);
%theta=acos(dotprod(sun_b,r./kron([1 1 1],dpmag)));

% Determine eclipse times
ec=find((rhop-rhos > theta));
out=ones(length(r),1);
out(ec)=zeros(length(ec),1);
