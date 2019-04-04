% ========================================================= %
%
% this script downloads satellite data over a given base
% receiver at a given time for a given length of time
% and visualizes the positions over time and the signal
% to noise ratio as a function of elevation
%
%
% The program will take roughly 45 seconds to run if 
% it downloads the data or roughly 15 seconds if the 
% data is current in the workspace. (Assuming its run 
% for a two hour interval)
%
%
% ========================================================= %
% 
% 10-28-18
% Author: Andrew Malatak
%
% References: 
% Dr. Moriba Jah Lecture notes, ASE 372N Fa18
%
% ========================================================= %





tic
gpsWeek = 1712;
gpsSec = 143800;
elMask = 15;
rRx = getAntLoc('WRW0');
cdtRx = 0;
tiflags = [0 0];



% retrieves satellite data for a given time and plots it
[satData, ionodata] = retrieveNavigationData(gpsWeek, gpsSec, 0);
[svIdVec, plotMat] = satmap(satData, rRx, elMask, gpsWeek, gpsSec, 0);

SVs = [];
for i = 1:length(plotMat(:,1))
	if plotMat(i, 3) > elMask*pi/180
		SVs(end+1) = svIdVec(i);
	end
end

for i = 1:length(SVs)
	[rhobar(i), Hrho(i, :)] = satpr(gpsWeek, gpsSec, cdtRx, rRx,satData(SVs(i)), ionodata, tiflags);
end


toc
