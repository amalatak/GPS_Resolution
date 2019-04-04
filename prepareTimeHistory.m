function [tRVec,obsValidMat,svIdVec,prMat,fDMat,thetaMat] = ...
      prepareTimeHistory(channelMat,epochStride)
% prepareTimeHistory : Generate a time history of measurements at each epoch
%                      and generate unified time vectors.
%
%
% INPUTS
%
% channelMat ---------- Matrix of data loaded from channel.mat
%
% epochStride --------- Stride length for data that will participate in the
%                       navigation solution.  Set to 1 to use all the data;
%                       set to 2 to skip every other epoch, etc.
%
%
% OUTPUTS
% 
% tRVec --------------- Structure of two Nt-by-1 vectors that jointly define
%                       the unique receiver time instants corresponding to all
%                       received observables; tRVec.w contains the GPS week
%                       and tRVec.s contains the GPS seconds of week.
%
% obsValidMat --------- Nt-by-Nsv matrix whose (ii,jj)th element indicates
%                       whether (1) or not (0) observables were valid for SVID
%                       svIdVec(jj) at time tRVec.s(ii).
%
% svIdVec ------------- Nsv-by-1 vector of unique SVIDs corresponding to
%                       SVs that were tracked at some point during the
%                       measurement interval tRVec.s(1) to tRVec.s(Nt).
%
% prMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
%                       measured pseudorange value for SVID svIdVec(jj) at
%                       time tRVec.s(ii), in meters.
%
% fDMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
%                       measured Doppler value for SVID svIdVec(jj) at time
%                       tRVec.s(ii), in Hz.
%
% thetaMat ------------ Nt-by-Nsv matrix whose (ii,jj)th element is the
%                       measured beat carrier phase value for SVID svIdVec(jj)
%                       at time tRVec.s(ii), in cycles.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys
%+==============================================================================+
  
columnDefinitions;
navConstants;
signalType = 0;

iidum = find(channelMat(:,signalTypeCol) == signalType);
channelMat = channelMat(iidum,:);
iidum = find(channelMat(:,prValidFlagCol) & ...
             ~channelMat(:,phaseErrorFlagCol));
tRVec.w = channelMat(iidum,ortWkCol);
tRVec.s = channelMat(iidum,ortWholeSecCol) + channelMat(iidum,ortFracSecCol);
[tRVec.s,I,J] = unique(tRVec.s);
tRVec.w = tRVec.w(I);
if(isempty(tRVec.s))
  error('No valid pseudoranges');
end
svIdVec = unique(channelMat(:,svIdCol));
Nsv = length(svIdVec);

% Isolate the valid observables corresponding to the times in tRVec
Nt = length(tRVec.s);
obsValidMat = zeros(Nt,Nsv);
prMat = zeros(Nt,Nsv);
fDMat = zeros(Nt,Nsv);
thetaMat = zeros(Nt,Nsv);
for jj=1:Nsv
  svId = svIdVec(jj);
  iidum = find(channelMat(:,svIdCol)==svId & ...
               channelMat(:,prValidFlagCol)==1 & ...
               ~channelMat(:,phaseErrorFlagCol));
  if(~isempty(iidum))
    sVMat = channelMat(iidum,:);
    tS.w = sVMat(:,ortWkCol);
    tS.s = sVMat(:,ortWholeSecCol) + sVMat(:,ortFracSecCol);
    [C,IA,IB] = intersect((tRVec.w - tRVec.w(1))*sec_in_week + tRVec.s, ...
                          (tS.w - tS.w(1))*sec_in_week + tS.s);
    obsValidMat(IA,jj) = ones(length(tS.s),1);
    prMat(IA,jj) = sVMat(IB,prCol);
    fDMat(IA,jj) = sVMat(IB,fpllCol);
    thetaMat(IA,jj) = sVMat(IB,thetahatCol);
  end
end
  

% Select indices separated by epochStride
iidum = 1:epochStride:length(tRVec.s);
tRVec.w = tRVec.w(iidum);
tRVec.s = tRVec.s(iidum);
obsValidMat = obsValidMat(iidum,:);
prMat = prMat(iidum,:);
fDMat = fDMat(iidum,:);
thetaMat = thetaMat(iidum,:);