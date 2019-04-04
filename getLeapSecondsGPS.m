function [leapSeconds] = getLeapSecondsGPS(gpsWeek,gpsSec)
% getLeapSecondsGPS : Returns the number of seconds by which GPS lead UTC at
%                     the specified input GPS epoch.
%
%
% INPUTS
%
% gpsWeek -------- The unambiguous GPS week number where zero corresponds
%                  to midnight on the evening of 5 January/morning of 6
%                  January, 1980.  By unambiguous is meant the full week
%                  count since the 1980 reference time (no rollover at
%                  1024 weeks).
%
% gpsSec --------- The GPS time of week expressed as GPS seconds from
%                  midnight on Saturday.  
%
%
% OUTPUTS
% 
% leapSeconds----- The number of seconds by which GPS time lead UTC at the
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
  
spw = 604800;  
gps = gpsWeek + gpsSec/spw;

% Find GPS times corresponding to the transition instant.  
 n = datenum(1997,6,30,0,0,0) + 1;
[gpsWeek12, gpsSec12] = utc2gps(n);
gps12 = gpsWeek12 + gpsSec12/spw; 
n = datenum(1998,12,31,0,0,0) + 1;
[gpsWeek13, gpsSec13] = utc2gps(n);  
gps13 = gpsWeek13 + gpsSec13/spw; 
n = datenum(2005,12,31,0,0,0) + 1;
[gpsWeek14, gpsSec14] = utc2gps(n);
gps14 = gpsWeek14 + gpsSec14/spw; 
 n = datenum(2008,12,31,0,0,0) + 1;
[gpsWeek15, gpsSec15] = utc2gps(n);
gps15 = gpsWeek15 + gpsSec15/spw; 
n = datenum(2012,6,30,0,0,0) + 1;
[gpsWeek16, gpsSec16] = utc2gps(n);
gps16 = gpsWeek16 + gpsSec16/spw;
n = datenum(2015,6,30,0,0,0) + 1;
[gpsWeek17, gpsSec17] = utc2gps(n);
gps17 = gpsWeek17 + gpsSec17/spw;
n = datenum(2016,12,31,0,0,0) + 1;
[gpsWeek18, gpsSec18] = utc2gps(n);
gps18 = gpsWeek18 + gpsSec18/spw;

leapSeconds = 0;
if gps < gps12
  leapSeconds = 11;
elseif (gps >= gps12 & gps < gps13)
  leapSeconds = 12;
elseif (gps >= gps13 & gps < gps14)
  leapSeconds = 13;
elseif (gps >= gps14 & gps < gps15)
  leapSeconds = 14; 
elseif (gps >= gps15 & gps < gps16)
  leapSeconds = 15;
elseif (gps >= gps16 & gps < gps17)
  leapSeconds = 16;
elseif (gps >= gps17 & gps < gps18)
  leapSeconds = 17;
else
  leapSeconds = 18;
end


  

