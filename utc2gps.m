function [gpsWeek, gpsSec] = utc2gps(n)
% utc2gps : Convert UTC time to GPS time expressed in GPS week and GPS
%           second of week.
%
%
% INPUTS
%
% n -------------- The UTC time and date expressed as a Matlab datenum. Use
%                  the Matlab function datenum() to generate n; use datestr()
%                  to render n in a standard format.
%
%
% OUTPUTS
%
% gpsWeek -------- The unambiguous GPS week number where zero corresponds to
%                  midnight on the evening of 5 January/morning of 6 January,
%                  1980. By unambiguous is meant the full week count since
%                  the 1980 reference time (no rollover at 1024 weeks).
%
% gpsSec --------- The GPS time of week expressed as GPS seconds from midnight
%                  on Saturday.
%
%+------------------------------------------------------------------------------+
% References:
%
% http://leapsecond.com/java/gpsclock.htm
% https://www.gps.gov/technical/icwg/IS-GPS-200J.pdf
%
% Author: Andrew Malatak
%+==============================================================================+

epoch = '01.06.1980';     		% epoch for gps time
format = 'mm.dd.yyyy';

datenumepc = datenum(epoch, format);    % formats epoch time

dd_snc_epc = n - datenumepc;            % days since epoch
sec_snc_epc = dd_snc_epc*86400;         % seconds since epoch

gpsWeek = floor(sec_snc_epc/604800);        % weeks since epoch
Secnoleap = sec_snc_epc - gpsWeek*604800;   % seconds since week start, no leap seconds

gpsSec = Secnoleap + getLeapSecondsUTC(n);  % seconds since start of week


%======================================================================%
% adds current and historical leap seconds - current as of (Oct. 2018) %



