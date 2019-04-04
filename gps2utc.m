function [n] = gps2utc(gpsWeek,gpsSec)
% gps2utc : Convert GPS time expressed in GPS week and GPS second of week to
%           UTC.
%
%
% INPUTS
%
% gpsWeek -------- The unambiguous GPS week number where zero corresponds to
%                  midnight on the evening of 5 January/morning of 6 January,
%                  1980. By unambiguous is meant the full week count since
%                  the 1980 reference time (no rollover at 1024 weeks).
%
% gpsSec --------- The GPS time of week expressed as GPS seconds from midnight
%                  on Saturday.
%
%
% OUTPUTS
%
% n -------------- The UTC time and date expressed as a Matlab datenum. Use
%                  the Matlab function datestr() to read n in a standard
%                  format.
%
%+------------------------------------------------------------------------------+
% References:
%
% http://leapsecond.com/java/gpsclock.htm
% https://www.gps.gov/technical/icwg/IS-GPS-200J.pdf
%
%
% Author: Andrew Malatak
%+==============================================================================+

lpsec = getLeapSecondsGPS(gpsWeek, gpsSec);  % number of leapseconds for the input gps time

gpssecnoleap = gpsSec - lpsec;               % gps seconds without leap seconds

secnoleap = gpssecnoleap + gpsWeek*604800;   % seconds since gps time epoch

ddsncepc = secnoleap/86400;                  % days since gps time epoch

epoch = '01.06.1980';                        % gps epoch
format = 'mm.dd.yyyy';
datenumepc = datenum(epoch, format);         % gps epoch as a datenum

n = datenumepc + ddsncepc;                   % utc datenum


%==================================================%
% getLeapSecondsGPS.m is current as of (Oct. 2018) %




