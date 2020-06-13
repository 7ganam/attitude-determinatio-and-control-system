function sat=mytle(tlefile, satname, mu)
%TLE2ORB  Get Keplerian elements from Two-Line Element data.
%   SAT = TLE2ORB([TLEFILE],satnema, mu), where SAT is an array of struct
%
%      SAT(1:N).oe
%      SAT(1:N).name
%      SAT(1:N).epoch
%
%   as described below. N is the number of satellites.
%
%   NAME is the satellite name,
%   EPOCH is in YMDhms format and
%   OE is the Keplerian orbital elements.
%   
%   Mean motion and it's derivatives are in rad/s, rad/s^2 and rad/s^3
%   Period is in seconds
%
%   The orbital elements have the format:
%      a [units from mu]   : semi-major axis
%      e []    : eccentricity
%      i [deg] : inclination
%      W [deg] : longitude of the ascending node
%      w [deg] : argument of perigee
%      M [deg] : mean anomaly at EPOCH
%      P [s]   : orbit period
% 
%   The data is extracted from the text file TLEFILE with 
%   Two-Line Element format of the satellite orbits.
%   To download current satellite data go to:
%
%      http://celestrak.com/NORAD/elements/
%
%   Then put the file(s) in the .../sat/ directory.
%
% Revision 2002-12-09, 2002-12-30, 2002-12-31,
%          2003-04-03, 2003-04-06, 2003-04-17,
%          2003-06-19.

fp=fopen(tlefile,'r');
if(fp == -1)
    error('Cannot open file');
end

%mu=my_const('G')*my_const('Me');
Pes=24*3600; % revolutions per solar day


%% scan file until satname found
while ~feof(fp)
    line=fgetl(fp);
    found = strcmp(line,satname);
    if found
        break
    end   
end
if ~found
    error('Cannot find that sat (may need to remove extra space after s/c name in file)')
end


%% read two lines
line1=fgetl(fp); 
line2=fgetl(fp); 
fclose(fp);

%% Information from LINE0
sat.name=deblank(satname);

%% Information from LINE1
sat.satid = line1(3:7);
sat.launchYear = str2double(line1(10:11));
sat.launchNumber = str2double(line1(12:14));
sat.launchpiece  = line1(15:17);
year=line1(19:20);
sat.epochDayOfYear=str2double(line1(21:32));
if str2double(year(1))>4
      syear=['19' year];
else
      syear=['20' year];
end
%dayToMin = 1440;
sat.epochYear=str2double(syear);
ep=datevec(datenum(sat.epochYear,0,sat.epochDayOfYear));
sat.epoch=ep;
sat.jdepoch = julian(ep(1),ep(2),ep(3),ep(4),ep(5),ep(6));
%sat.ballisticCoeff   = str2double(line1(34:43));
sat.n0               = str2double(line2(53:63))*2*pi/Pes;
sat.n0Dot            = str2double(line1(34:43))*2*2*pi/(Pes^2);
sat.n0DDot           = Str2NumE(line1(45:52))*6*2*pi/(Pes^3)/1e5;
sat.bStar            = Str2NumE(line1(54:61))/1e5;
%sat.radPressure      = Str2NumE(line1(54:61))/1e5;

%% Information from LINE2
sat.rpd=str2double(line2(53:63));       % revs per day
P=Pes/sat.rpd;                          % [sec] orbital period 
OE(1)=(mu/(sat.rpd*2*pi/Pes)^2)^(1/3);  % semimajor axis
OE(2)=str2double(['.' line2(27:33)]);   % eccentricity
OE(3)=str2double(line2(9:16));          % [deg] inclination
OE(4)=str2double(line2(18:25));         % [deg] longitude (or right ascension) of the ascending node 
OE(5)=str2double(line2(35:42));         % [deg] argument of periapsis (or perigee)
OE(6)=str2double(line2(44:51));         % [deg] mean anomaly at epoch
OE(7)=P;                                % [sec] orbit period

sat.oe=OE;
sat.p=P;

function x = Str2NumE( s )
x  = str2double(s(1:6));
xE = str2double(s(7:8));
x  = x*10^xE;

