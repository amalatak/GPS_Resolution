function [SVIDmat, plotMat] = satmap(satdata,rRxEcef,elMaskDeg,gpsWeek,gpsSecVec,plotFlag)
% satmap : Generate plotting data for the SVs above a particular receiver
% position over a span of time.
%
%
% INPUTS
%
% satdata ------ Ephemeris structure array; see getephem header for details.
%
% rRxEcef ------ 3-by-1 receiver ECEF position, in meters.
%
% elMaskDeg ---- Elevation mask angle, in degrees.
%
% gpsWeek ------ GPS week number corresponding to the first element in
% gpsSecVec.
%
% gpsSecVec ---- Nt-by-1 vector of GPS seconds of week over which an SV trace
% is desired. Entries exceeding 7*86400 indicate that the time
% interval spans a GPS week boundary.
%
% plotFlag ----- Indicates whether (1) or not (0) to generate a sky plot of
% SVs.
%
%
% OUTPUTS
%
% svIdVec ------ Nsv-by-1 vector of unique SV identification numbers for the
% SVs that were above the elevation mask angle at any time
% during the interval.
%
% plotMat ------ Nt*Nsv-by-4 matrix of data corresponding to SVs in svIdVec,
% arranged as
%
% [svID gpsSec elRad azRad;
% svID gpsSec elRad azRad;
% svID gpsSec elRad azRad]
%
% where svId is the SV identification number, gpsSec is GPS
% time in seconds of week, and elRad and azRad are the SV
% elevation and azimuth angles, in radians. plotMat is
% composed of Nt batches of Nsv rows each, where the ith batch
% corresponds to time gpsSecVec(i). Each batch has data for
% all Nsv svIDs in svIdVec. Note that some elRad values may be
% below the elevation mask angle: this simply means that the
% corresponding svID rose or set during the data interval.
%
%+------------------------------------------------------------------------------+
% References:
% Dr. Moriba Jah Lecture notes, ASE 372N Fa18
%
% Author: Andrew Malatak
%+==============================================================================+

% variable initialization

plotmatsize = length(satdata)*length(gpsSecVec);
plotMat = zeros(plotmatsize, 4);

SVIDsize = length(gpsSecVec);
SVIDmat = zeros(SVIDsize, 1);

emptyrows = zeros(plotmatsize, 1);

count = 0;


% Obtains Satellite plotMat matrix and SVID matrix
% after checking for good data. 
for j = 1:length(gpsSecVec)
    for i = 1:length(satdata)
            count = count + 1;
	if ~isempty(satdata(i).SVID)
	    recef = satloc(gpsWeek, gpsSecVec(j), satdata(i));
	    [el, az] = satelaz(recef, rRxEcef);
	    svId = satdata(i).SVID;
            plotMat(count, 1:4) = [svId, gpsSecVec(j), el, az];
            SVIDmat(j, i) = svId;
        else 
            emptyrows(count) = 1;
	end
    end
end

% eliminate empty satellite rows
for i = 1:length(emptyrows)
    if emptyrows(i)
        plotMat(i, 1:4) = [];
    end
end



% if plot flag is true, plot the satellites above the receiver location
if plotFlag
    plotsat(plotMat, gpsWeek, elMaskDeg*pi()/180);
end
