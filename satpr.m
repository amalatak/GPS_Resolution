function [rhoBar, Hrho] = ...
	satpr(gpsWeekRx, gpsSecRx, cdtRx, rRxEcef, sd, ionodata, tiFlags)
% satpr : Calculate the modeled pseudorange between the satellite specified in
% sd and the receiver located at rRxEcef at the given time of signal
% reception.
%
%
% INPUTS
%
% gpsWeekRx -- GPS week of signal reception event in receiver time.
%
% gpsSecRx --- GPS seconds of week of signal reception event in receiver
% time.
%
% cdtRx ------ Receiver clock error scaled by the speed of light: cdtRx =
% c*dtRx. True time t is related to receiver time tRx by t = tRx
% - dtRx.
%
% rRxEcef ---- 3-by-1 ECEF coordinates of the receiver antenna’s phase center,
% in meters.
%
% sd --------- Ephemeris structure array for a single SV. Let ii be the
% numerical identifier (PRN identifier) for the SV whose location
% is sought. Then sd = satdata(ii) has the following fields:
%
% SVID - satellite number
% health - satellite health flag (0 = healthy; otherwise unhealthy)
% we - week of ephemeris epoch (GPS week, unambiguous)
% te - time of ephemeris epoch (GPS seconds of week)
% wc - week of clock epoch (GPS week)
% tc - time of clock epoch (GPS seconds of week)
% e - eccentricity (unitless)
% sqrta - sqrt of orbit semi-major axis (m^1/2)
% omega0 - argument of perigee (rad.)
% M0 - mean anomaly at epoch (rad.)
% L0 - longitude of ascending node at beginning of week (rad.)
% i0 - inclination angle at epoch (rad.)
% dOdt - longitude rate (rad / sec.)
% dn - mean motion difference (rad / sec.)
% didt - inclination rate (rad / sec.)
% Cuc - cosine correction to argument of perigee (rad.)
% Cus - sine correction to argument of perigee (rad.)
% Crc - cosine correction to orbital radius (m)
% Crs - sine correction to orbital radius (m)
% Cic - cosine correction to inclination (rad.)
% Cis - sine correction to inclination (rad.)
% af0 - 0th order satellite clock correction (s)
% af1 - 1st order satellite clock correction (s / s)
% af2 - 2nd order satellite clock correction (s / s^2)
% TGD - group delay time for the satellite (s)
%
% ionodata -- Ionospheric data structure array with the following fields:
%
% alpha0, alpha1, alpha2, alpha3 - power series expansion coefficients
% for amplitude of ionospheric TEC
%
% beta0, beta1, beta2, beta3 - power series expansion coefficients for
% period of ionospheric plasma density cycle
%
% tiFlags --- 2-by-1 vector of flags that indicate whether to enable
% tropospheric delay correction (tiFlag(1)) or ionospheric
% delay correction (tiFlags(2)).
%
%
% OUTPUTS
%
% rhoBar ---- Modeled pseudorange between the receiver at time of the signal
% reception event (the time specified by gpsWeekRx and gpsSecRx)
% and the satellite at time of the signal transmission event, in
% meters.
%
% Hrho ------ 1 x 4 Jacobian matrix containing partials of pseudorange with
% respect to rRxEcef and c*dtRx
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:
%+==============================================================================+

c = 299792458;
GM = 3.986005e14;

tr = gpsWeekRx*7*24*3600 + gpsSecRx;



E = fzero(@(x) x - sd.e*sin(x) - sd.M0, sd.M0);

tc = sd.wc*7*24*3600 + sd.tc;
tgd = sd.TGD;
dtecc = -2*sd.e*sin(E)*sqrt(GM)*sd.sqrta/(c^2);

dtubar = cdtRx/c;

% Towlt

Towlt = .075;
drst = inf;
rstold = 0;

while drst > .00001
	tbar = tr - dtubar - Towlt;
	tbarwk = floor(tbar/(7*24*3600));
	tbarsc = tbar - tbarwk*7*24*3600;
	rst = satloc(tbarwk, tbarsc, sd);
	Towlt = (1/c)*norm(rRxEcef - rst);
	drst = rst - rstold;
	rstold = rst;
end

Towlt*c;

tbar = tr - dtubar - Towlt;

dts = sd.af0 + sd.af1*(tbar - tc) + sd.af2*(tbar - tc)^2 + dtecc - tgd;

OmegaE = 7.2921151467e-5;

rstecef = R3(OmegaE*Towlt)*rst;

rhoBar = norm(rRxEcef - rstecef) + cdtRx - c*dts;


if size(rRxEcef) == size(rstecef)
	Hi = (rRxEcef - rstecef)/norm(rRxEcef - rstecef);
else
	Hi = (rRxEcef' - rstecef)/norm(rRxEcef' - rstecef);
end
Hrho = [Hi' 1];

cf = 1575.42e6; % carrier frequency
IonoModel = 'broadcast'; % Klobuchar Model
tGPS.seconds = sd.tc; tGPS.week = sd.wc;

delTauG = getIonoDelay(ionodata, cf, rRxEcef, rst, tGPS, IonoModel);
IonoDelay = delTauG*c; % Ionospheric Delay

rst; % off a little ===========================================================

tropomodel = 'Saastamoinen_MSP_Neill';
delTauNa = getTropoDelay(rRxEcef, rst, tGPS, tropomodel);
TropoDelay = delTauNa*c; % Tropospheric Delay


rhoBar = rhoBar + IonoDelay + TropoDelay;
