function [solutionMat,tRVecSolution,tTrueVecSolution,...
residualVec,badSvIdVec,thetaNominalMat] = ...
performNavigationSolution(tRVec,svIdVec,obsValidMat,prMat,fDMat,thetaMat,...
iiStart,iiStop,rRxAppx, ...
satdata,ionodata,svIdAllow,lambda,tiFlags)

% performNavigationSolution : Perform the navigation solution using nonlinear
% least squares techniques.
%
%
% INPUTS
%
% tRVec --------------- Structure of two Nt-by-1 vectors that jointly define
% the unique receiver time instants corresponding to all
% received observables; tRVec.w contains the GPS week
% and tRVec.s contains the GPS seconds of week.
%
% svIdVec ------------- Nsv-by-1 vector of unique SVIDs corresponding to
% SVs that were tracked at some point during the
% measurement interval spanned by tRVec.
%
% obsValidMat --------- Nt-by-Nsv matrix whose (ii,jj)th element indicates
% whether (1) or not (0) observables were valid for SVID
% svIdVec(jj) at time tRVec.s(ii).
%
% prMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured pseudorange value for SVID svIdVec(jj) at
% time tRVec.s(ii), in meters.
%
% fDMat --------------- Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured Doppler value for SVID svIdVec(jj) at time
% tRVec.s(ii), in Hz.
%
% thetaMat ------------ Nt-by-Nsv matrix whose (ii,jj)th element is the
% measured beat carrier phase value for SVID svIdVec(jj)
% at time tRVec.s(ii), in cycles.
%
% iiStart,iiStop ------ Start and stop indices into tRVec between which a
% navigation solution will be computed. These indices
% are used to isolate a specific section of data. The
% number of solutions Ns will nominally be Ns = iiStop -
% iiStart + 1 but will be less if for some epochs a
% solution could not be computed (e.g., too few SVs).
%
% rRxAppx ------------- 3-by-1 approximate location of receiver antenna, in
% ECEF meters.
%
% satdata, ionodata --- SV and ionospheric navigation data parameters.
% See getephem.m for details.
%
% svIdAllow ----------- Vector of SVIDs for SVs that are allowed to
% participate in the navigation solution.
%
% lambda -------------- Nominal wavelength of GNSS signal, in meters.
%
%
% tiFlags ------------- 2-by-1 vector of flags that indicate whether to enable
% tropospheric delay correction (tiFlag(1)) or
% ionospheric delay correction (tiFlags(2)).
%
%
% OUTPUTS
%
% solutionMat --------- Ns-by-4 matrix of navigation solutions whose rows
% are formatted as [xEcef,yEcef,zEcef,c*deltR].
%
% tRVecSolution ------- Structure of two Ns-by-1 vectors that jointly define
% the unique receiver time instants corresponding to
% navigation solutions in solutionMat; tRVecSolution.w
% contains the GPS week and tRVecSolution.s contains the
% GPS seconds of week.
%
% tTrueVecSolution ---- Structure of two Ns-by-1 vectors that jointly define
% the unique receiver time instants corresponding to
% navigation solutions in solutionMat. This vector is
% an estimate of ’true’ GPS time for each
% solution. tTrueVecSolution.w contains the GPS week and
% tTrueVecSolution.s contains the GPS seconds of week.
%
% residualVec --------- Ns-by-1 vector of maximum least-squares residuals
% magnitudes from each solution epoch, in meters.
%
% badSvIdVec ---------- Ns-by-1 vector of SVIDs corresponding to the SV with
% the worst residual at each solution epoch.
%
% thetaNominalMat ----- Ns-by-Nsv vector of nominal beat carrier phase values,
% in the same arrangement as those in thetaMat, based
% solely on the SV-to-rRxAppx range evaluated at
% receiver time, in cycles.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:
%+==============================================================================+

c = 299792458;
sigma = 2;
xnew = [rRxAppx; 0];

dx = inf;
dp = 0;

convflag = inf;
initsol = [0; 0; 0; 0];
flag = false;

for time = iiStart:iiStop

	fprintf('Solving for time %d\n\n', time)

	timeindex = time - iiStart + 1;

	while norm(convflag) > 1000 || flag
		flag = false;
		while (dx > 1e-6)


			xbar = xnew;
			zprime(timeindex,:) = prMat(time, :);

			gpsWeek = tRVec.w(time); gpsSec = tRVec.s(time);
	                [sd, ionodata] = retrieveNavigationData(gpsWeek, gpsSec, 0);


			for SV = 1:length(prMat(1,:))

				if obsValidMat(time, SV)

					SVID = svIdVec(SV);
					gpsWeek = tRVec.w(time); gpsSec = tRVec.s(time);
					[p(SV,1), H(SV,1:4)] = satpr(gpsWeek, gpsSec, xbar(4), ...
						 xbar(1:3), satdata(SVID), ionodata, [0 0]);
					pr(SV, 1) = prMat(time, SV);
					dp(SV) = pr(SV) - p(SV);

					if dp(SV) >= max(dp)
						residualVec = dp(SV);
						badSvIdVec = SVID;
					end
				end
			end

			Hstar = (1/sigma)*H;
			dp = pr - p;
			dpstar = (1/sigma)*dp;
			dx = inv(Hstar'*Hstar)*Hstar'*dpstar;
			xnew = xbar + dx;
		end

		[lat, lon, h] = ecef2lla(xnew);

		dx = inf;
		solutionMat(timeindex, :) = xnew;
		tRVecSolution(timeindex).w = gpsWeek
		tRVecSolution(timeindex).s = gpsSec

		tTrueVecSolution(timeindex).w = gpsWeek
		tTrueVecSolution(timeindex).s = gpsSec + xnew(4)/c

		convflag = xnew - initsol;
		initsol = xnew;
	end
	flag = true;
end


thetaNominalMat = [];
end

