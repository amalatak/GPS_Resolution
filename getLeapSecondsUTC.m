function [leapSeconds] = getLeapSecondsUTC(n)
% getLeapSecondsUTC : Returns the number of seconds by which GPS lead UTC at
%                     the specified input UTC epoch.
%
%
% INPUTS
%
% n -------------- The UTC time and date expressed as a Matlab datenum.  Use
%                  the Matlab function datestr() to read n in a standard
%                  format.
%
%
% OUTPUTS
% 
% leapSeconds ---- The number of seconds by which GPS time lead UTC at the
%                  specified epoch.    
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys
%+==============================================================================+
  
  
% Notes: These are the leap seconds that came into effect after midnight
% UTC on the indicated days.
%
% June     30, 1997     	12 seconds
% December 31, 1998		13 seconds
% December 31, 2005		14 seconds
% December 31, 2008		15 seconds
% June     30, 2012     	16 seconds
% June     30, 2015             17 seconds
% December 31, 2016             18 seconds 


% Find UTC times corresponding to the transition instant.
n12 = datenum(1997, 6,30,0,0,0) + 1;
n13 = datenum(1998,12,31,0,0,0) + 1;
n14 = datenum(2005,12,31,0,0,0) + 1;
n15 = datenum(2008,12,31,0,0,0) + 1;
n16 = datenum(2012, 6,30,0,0,0) + 1;
n17 = datenum(2015, 6,30,0,0,0) + 1;
n18 = datenum(2016,12,31,0,0,0) + 1;

leapSeconds = 0;
if n < n12
  leapSeconds = 11;
elseif (n >= n12 & n < n13)
  leapSeconds = 12;
elseif (n >= n13 & n < n14)
  leapSeconds = 13;
elseif (n >= n14 & n < n15)
  leapSeconds = 14;
elseif (n >= n15 & n < n16)
  leapSeconds = 15;
elseif (n >= n16 & n < n17)
  leapSeconds = 16;
elseif (n >= n17 & n < n18)
  leapSeconds = 17;
else
  leapSeconds = 18;
end




