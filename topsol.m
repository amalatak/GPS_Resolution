% topsol.m
%
% Top-level script for computing a navigation solution from
% channel.mat-formatted observables

clear;clc;fclose('all');
%----- Setup
% Elevation mask angle, in degrees
elMaskDeg = 0;
% Exclude the following svIds from participating in the solution. SVs marked
% unhealthy will be automatically added to this list.
svIdExclude = [];
% Set start and stop indices
iiStart = 1; iiStop = 288;
% Set approximate receiver antenna position
[ra,va,ya] = getAntLoc('RHO1');
%rRAppx = ra + va*(2016.8 - ya);
rRAppx = getAntLoc('Fiction_Island')
% Enable/disable tropospheric and ionospheric corrections. 
tiFlags(1) = 1; tiFlags(2) = 1;
% Name of output file for Google Earth rendering
kmlFileOutName = 'kmlOut.kml';
% Change to -3 if data are so recent that no matching ephemeris is found
ephemHourOffset = 0;
% Ephemeris/iono data and elevation masking refresh interval
subSolutionBlockIntervalSec = 50;
% Observation epoch stride length
epochStride = 1;

%----- Load and prepare data
navConstants; load(['channel.mat']); channelMat = channel';
[tRVec,obsValidMat,svIdVec,prMat,fDMat,thetaMat] = ...
    prepareTimeHistory(channel,epochStride);

delii = ceil(subSolutionBlockIntervalSec/min(diff(tRVec.s)));

%----- Iterate on sub-solution blocks
solutionMat = []; tRVecSolution.w = []; tTrueVecSolution.w = [];
tRVecSolution.s = []; tTrueVecSolution.s = [];
residualVec = []; badSvIdVec = []; thetaNominalMat = [];
iiStop = min(iiStop,length(tRVec.s));
nii = length(iiStart:delii:iiStop);
fprintf('Processing ... \n');
for(mm = 1:nii) fprintf('='); end; fprintf('\n');
for iiA = iiStart:delii:iiStop
  % Refresh ephemeris/iono data and masking
  iiB = min(iiStop,iiA + delii - 1); iiM = round(mean([iiA,iiB]));
  [satdata, ionodata] = retrieveNavigationData(tRVec.w(iiM),tRVec.s(iiM),...
                                               ephemHourOffset,'navFiles');
  [svIdAllow,] = satmap(satdata,rRAppx,elMaskDeg,tRVec.w(iiM),tRVec.s(iiM),1);

  % Exclude unhealthy PRNs
  svIdExcludeLocal = svIdExclude;
  for ii=1:length(satdata)
    if(satdata(ii).health) svIdExcludeLocal = [svIdExcludeLocal,ii]; end
  end

  svIdAllow = setdiff(svIdAllow,svIdExcludeLocal(:)); fprintf('=');
  % Perform nav solution

  [solutionMatD,tRVecSolutionD,tTrueVecSolutionD,...
   residualVecD,badSvIdVecD,thetaNominalMatD] = ...
      performNavigationSolution(tRVec,svIdVec,obsValidMat,prMat,fDMat,...
                                thetaMat,iiA,iiB,rRAppx,satdata,ionodata,...
                                svIdAllow,lambdaGPSL1,tiFlags);

  % Store data
  solutionMat = [solutionMat;solutionMatD];
  tRVecSolution.w = [tRVecSolution.w;tRVecSolutionD.w];
  tRVecSolution.s = [tRVecSolution.s;tRVecSolutionD.s];
  tTrueVecSolution.w = [tTrueVecSolution.w;tTrueVecSolutionD.w];
  tTrueVecSolution.s = [tTrueVecSolution.s;tTrueVecSolutionD.s];
  residualVec = [residualVec;residualVecD];
  badSvIdVec = [badSvIdVec;badSvIdVecD];
  thetaNominalMat = [thetaNominalMat;thetaNominalMatD];
end

%----- Solve for average position
rRSolAvg = mean(solutionMat(:,1:3),1)';
[latAvg,lonAvg,altAvg] = ecef2lla(rRSolAvg);
[Recef_to_enu]=ecef2enu(latAvg,lonAvg);
dSolution = [solutionMat(:,1)-rRSolAvg(1),solutionMat(:,2)-rRSolAvg(2),...
             solutionMat(:,3)-rRSolAvg(3)];
dSolutionENU = (Recef_to_enu*dSolution')';
dSolutionAvgFromPriorENU = Recef_to_enu*(rRAppx - rRSolAvg);

%{
%----- Plot results
iiVec = iiStart:(iiStart + length(tTrueVecSolution.s) - 1);
figure(5);clf;
plot(dSolutionENU(:,1),dSolutionENU(:,2), '.');
hold on;
plot(dSolutionAvgFromPriorENU(1),dSolutionAvgFromPriorENU(2), 'r.', ...
     'markersize', 30);
xlabel('East (meters)'); ylabel('North (meters)');
title('Horizontal spread of navigation solution points with origin at the mean');
grid on; axis equal
figure(6);clf;
plot(iiVec,dSolutionENU(:,3));
hold on;
plot(iiVec,ones(length(iiVec),1)*dSolutionAvgFromPriorENU(3), 'r', ...
     'linewidth', 2);
xlabel('Solution Number'); ylabel('displacement (meters)');
title('Vertical displacement from mean');
figure(7);clf;
plot(iiVec,residualVec);
xlabel('Solution Number'); ylabel('max \Delta z (meters)');
title('Worst case residual at each solution epoch');
figure(8);clf;
plot(iiVec,badSvIdVec,'.');
title('SV with largest residual at each epoch');
xlabel('Solution Number');
grid on;

%}
%----- Print solution
fprintf('\n\nMean navigation solution:\n');
fprintf(' Latitude (deg)   : %+f\n', latAvg*180/pi);
fprintf(' Longitude (deg)  : %+f\n', lonAvg*180/pi);
fprintf(' Altitude (meters): %+f\n', altAvg);
fprintf('Mean epoch spacing  (seconds): %f\n', ...
        mean(diff((tTrueVecSolution.w - tTrueVecSolution.w(1))*sec_in_week + ...
                  tTrueVecSolution.s)));
fprintf('Total time interval (seconds): %f\n', ...
        (tTrueVecSolution.w(end) - tTrueVecSolution.w(1))*sec_in_week + ...
        tTrueVecSolution.s(end) - tTrueVecSolution.s(1));
fprintf('Solution distance from initial guess (meters): %f\n', ...
        norm(dSolutionAvgFromPriorENU));
fprintf('Solution horizontal distance from initial guess (meters): %f\n', ...
        norm(dSolutionAvgFromPriorENU(1:2)));
fprintf('Solution vertical distance from initial guess   (meters): %f\n', ...
        abs(dSolutionAvgFromPriorENU(3)));

% %----- Generate Google Earth KML file
% pMat = solutionMat(:,1:3);
% tVec.week = tTrueVecSolution.w;
% tVec.sec = tTrueVecSolution.s;
% genKmlFile(pMat,tVec,,kmlFileOutName);
% st = ['! cp ' kmlFileOutName ' /home/todd '];
% eval(st);

% %----- Create navsol.mat for use in CDGNSS processing; save in obsPath
% columnDefinitions;
% wkInvalid = 9999;
% nS = length(tRVecSolution.s);
% navsol = zeros(nS,12);
% navsol(:,ns_ortWkCol) = tRVecSolution.w;
% navsol(:,ns_ortWholeSecCol) = floor(tRVecSolution.s);
% navsol(:,ns_ortFracSecCol) = tRVecSolution.s - floor(tRVecSolution.s);
% navsol(:,4:7) = solutionMat;
% navsol(:,12) = 2*ones(nS,1);
% navsol = navsol';
% save([obsPath '/navsol.mat'], 'navsol');


