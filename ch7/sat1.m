function out=sat1(s,e);
%function out=sat1(s,e);
% This function outputs the saturation function.
% This function processes a vector of input x. 
% The output is defined as:
%      sat(x) = 1         if s > e
%      sat(x) = s/e       if abs(s) <= e
%      sat(x) = -1        if s < -e

if (s > e),
 out = 1;
end

if (abs(s) <= e),
 out = s/e;
end

if (s < -e),
 out = -1;
end