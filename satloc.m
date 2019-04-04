function [rSvEcef, vSvEcef] = satloc(gpsWeek, gpsSec, sd)
% satloc : Return satellite location and velocity expressed in and relative to
%	   the ECEF reference frame.
%
%
% INPUTS
%
% gpsWeek ---- Week of true time at which SV location and velocity are
%	       desired.
%
% gpsSec ----- Seconds of week of true time at which SV location and velocity
%	       are desired.
%
% sd --------- Ephemeris structure array for a single SV. Let ii be the
%	       numerical identifier (PRN identifier) for the SV whose location
%	       is sought. Then sd = satdata(ii). sd has the following
%	       fields:
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
%
% OUTPUTS
%
% rSvEcef ---- Location of SV at desired time expressed in the ECEF reference
%	       frame (m).
%
% vSvEcef ---- Velocity of SV at desired time relative to the ECEF reference
%	       frame and expressed in the ECEF reference frame (m/s). NOTE:
% 	       vSvEcef is NOT inertial velocity, e.g., for geostationary SVs
% 	       vSvEcef = 0.
%
%+------------------------------------------------------------------------------+
% References:
%
%
%
%
% Author: Andrew Malatak
%+==============================================================================+



%============================================%
          %%% Other Constants %%%

Wearth =  7.2921151467e-5;
GM = 3.986005e14;

%============================================%
      %%% Position Calculation %%%

A = sd.sqrta^2;              % semi major axis
n0 = sqrt(GM/(A^3));         % mean motion
n = n0 + sd.dn;              % corrected mean motion

tk = 604800*(gpsWeek - sd.we) + gpsSec - sd.te; % time from ephemeris reference epoch
Mk = sd.M0 + n*tk;                              % Mean anomaly
Ek = fzero(@(x) x - sd.e*sin(x) - Mk, Mk);      % Eccentric Anomaly

numer = sqrt(1 - sd.e^2)*sin(Ek);
denom = (cos(Ek) - sd.e);
vk = atan2(numer, denom);         %True Anomaly
phik = vk + sd.omega0;            %Argument of Latitude

duk = sd.Cus*sin(2*phik) + sd.Cuc*cos(2*phik);   % Second Harmonic Perturbations
drk = sd.Crs*sin(2*phik) + sd.Crc*cos(2*phik);
dik = sd.Cis*sin(2*phik) + sd.Cic*cos(2*phik);

uk = phik + duk;                         % Corrected Argument of Latitude
rk = A*(1 - sd.e*cos(Ek)) + drk;         % Corrected Radius
ik = sd.i0 + dik + sd.didt*tk;           % Corrected Inclination

xkp = rk*cos(uk);		         % Perifocal Position
ykp = rk*sin(uk);

Ok = sd.L0 + (sd.dOdt - Wearth)*tk - Wearth*sd.te;  % Corrected Longitude of Asencding Node


% ECEF Position %

xk = xkp*cos(Ok) - ykp*cos(ik)*sin(Ok);
yk = xkp*sin(Ok) + ykp*cos(ik)*cos(Ok);
zk = ykp*sin(ik);

rSvEcef = [xk; yk; zk];



%============================================%
      %%% Velocity Calculation %%%

%% all variables with 'dot' are derivatives of variables
%% in the position calculation section

Mkdot = n;
Ekdot = Mkdot/(1 - sd.e*cos(Ek));
vkdot = sin(Ek)*Ekdot*(1 + sd.e*cos(vk))/((1 - cos(Ek)*sd.e)*sin(vk));

phikdot = vkdot;

dukdot = 2*phikdot*(sd.Cus*cos(2*phik) - sd.Cuc*sin(2*phik));
drkdot = 2*phikdot*(sd.Crs*cos(2*phik) - sd.Crc*sin(2*phik));
dikdot = 2*phikdot*(sd.Cis*cos(2*phik) - sd.Cic*sin(2*phik));

ukdot = phikdot + dukdot;
rkdot = A*sd.e*sin(Ek)*Ekdot + drkdot;
ikdot = sd.didt + dikdot;

xkpdot = rkdot*cos(uk) - ykp*ukdot;
ykpdot = rkdot*sin(uk) + xkp*ukdot;

Okdot = sd.dOdt - Wearth;

% ECEF velocity %

xkdot = xkpdot*cos(Ok) - ykpdot*cos(ik)*sin(Ok) + ykp*sin(ik)*sin(Ok)*ikdot - yk*Okdot;
ykdot = xkpdot*sin(Ok) + ykpdot*cos(ik)*cos(Ok) - ykp*sin(ik)*ikdot*cos(Ok) + xk*Okdot;
zkdot = ykpdot*sin(ik) + ykp*cos(ik)*ikdot;


vSvEcef = [xkdot; ykdot; zkdot];


