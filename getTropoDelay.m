function [delTauNa] = getTropoDelay(rRx,rSv,tGPS,model)
% getTropoDelay : Return a model-based estimate of the neutral atmospheric
%                 delay experienced by a GNSS signal as it propagates from a
%                 GNSS SV to the antenna of a terrestrial GNSS receiver.
%
% INPUTS
%
% rRx ------------ A 3-by-1 vector representing the receiver antenna position
%                  at the time of receipt of the signal, expressed in meters
%                  in the ECEF reference frame.
%
% rSv ------------ A 3-by-1 vector representing the space vehicle antenna
%                  position at the time of transmission of the signal,
%                  expressed in meters in the ECEF reference frame.
%
% tGPS ----------- A structure containing the true GPS time of receipt of
%                  the signal.  The structure has the following fields:
%
%                  week -- unambiguous GPS week number
%
%                  seconds -- seconds (including fractional seconds) of the
%                  GPS week
%
% model ---------- A string identifying the model to be used in the
%                  computation of the ionospheric delay:
%
%                  Saastamoinen_MSP_Neill --- Saastamoinen model with
%                  modeled surface parameters (temperature, pressure) and
%                  the Neill mapping function.
%
%                  Other models TBD ...
%
% OUTPUTS
%
% delTauNa ------- Modeled scalar excess neutral atmospheric delay experienced
%                  by the GNSS signal, in seconds.
% 
%+------------------------------------------------------------------------------+
% Saastamoinen 1971, Davis 1986, Neill 1995.
%
%+==============================================================================+
% Copyright © 2006 Isaac Miller, all rights reserved.  No portion of this code
% may be reproduced, reused, or distributed in any form without prior written
% permission of the author.
  
if(strcmp(model,'Saastamoinen_MSP_Neill')~=1)
  error('Unrecognized model.');
end

load NMFcoeffs.mat;
navConstants;

% receiver geodetic latitude and longitude, in radians
[lat,lon,alt] = ecef2lla(rRx);
% satellite elevation and azimuth
[sEL, sAZ] = satelaz(rSv,rRx);

%set hydrostatic and wet zenith neutral delays (m)
%average hydrostatic zenith delay value (Neill 1995):
%znd_hydro = 2.3;
%average wet zenith delay value (Neill 1995):
%znd_wet = 0.1;
%values derived from thermodynamics:
%current surface pressure in Pascals
Po = 101325;
%current surface temperature To in Kelvins
To = 294.261111;
%current relative humidity RH (pct.)
RH = 60;
%hydrostatic zenith delay value derived from thermodynamics (Davis 1986):
%latitude correction factor to gravity (ellipsoidal)
lat_fact = (1 - 0.00266*cos(2*lat) - 0.00028*alt/1000);
znd_hydro = 0.000022768*Po / lat_fact;
%wet zenith delay value derived from thermodynamics (Saastamoinen 1971):
%Saastamoinen surface gravity approximation (m/s^2)
g = 9.784*lat_fact;
%empirical approximation to the saturation pressure of water (Pa)
pstarw = 3386.39*((0.00738*(To-273.15) + 0.8072)^8 ...
    - 0.000019*abs(1.8*(To-273.15) + 48) + 0.001316);
%partial pressure of water (Pa)
pw = RH/100 * pstarw;
%value derived from thermodynamics (Saastamoinen 1971):
znd_wet = (77.624e-8*scgRd/(4*g)*(1 - scgRd/scgRw) ...
    - scgRd*saast_cw/(4*g) ...
    + saast_cwp/(4*g/scgRd + stBeta)/To) * pw;
%slant factor based on geodetic latitude, altitude, and time (Neill 1995):
phs = 0;
if (lat < 0)
    %southern hemisphere is antisymmetric to northern
    phs = pi;
end
alat = abs(lat);
if (alat <= 0.261799387799149)
    %less than 15o latitude, use the 15o terms
    amp = NMFHamp(:, 1);
    avg = NMFHavg(:, 1);
    wet = NMFW(:, 1);
elseif (alat <= 0.523598775598299)
    %less than 30o, interpolate
    w = 1 - (alat - 0.261799387799149)/0.261799387799149;
    amp = w*NMFHamp(:, 1) + (1 - w)*NMFHamp(:, 2);
    avg = w*NMFHavg(:, 1) + (1 - w)*NMFHavg(:, 2);
    wet = w*NMFW(:, 1) + (1 - w)*NMFW(:, 2);
elseif (alat <= 0.785398163397448)
    %less than 45o, interpolate
    w = 1 - (alat - 0.523598775598299)/0.261799387799149;
    amp = w*NMFHamp(:, 2) + (1 - w)*NMFHamp(:, 3);
    avg = w*NMFHavg(:, 2) + (1 - w)*NMFHavg(:, 3);
    wet = w*NMFW(:, 2) + (1 - w)*NMFW(:, 3);
elseif (alat <= 1.0471975511966)
    %less than 60o, interpolate
    w = 1 - (alat - 0.785398163397448)/0.261799387799149;
    amp = w*NMFHamp(:, 3) + (1 - w)*NMFHamp(:, 4);
    avg = w*NMFHavg(:, 3) + (1 - w)*NMFHavg(:, 4);
    wet = w*NMFW(:, 3) + (1 - w)*NMFW(:, 4);
elseif (alat <= 1.30899693899575)
    %less than 75o, interpolate
    w = 1 - (alat - 1.0471975511966)/0.261799387799149;
    amp = w*NMFHamp(:, 4) + (1 - w)*NMFHamp(:, 5);
    avg = w*NMFHavg(:, 4) + (1 - w)*NMFHavg(:, 5);
    wet = w*NMFW(:, 4) + (1 - w)*NMFW(:, 5);
else
    %greater than 75o latitude, use the 75o terms
    amp = NMFHamp(:, 5);
    avg = NMFHavg(:, 5);
    wet = NMFW(:, 5);
end
ht = NMFHht;
%calculate day of year phase (time since DOY 28)
%DOY 28, 2006 is Jan. 28, week 1359, day 7
doyp = (tGPS.week - 1359)*7 + (tGPS.seconds - 518400)/86400;
doyp = 2*pi*doyp/365.25 + phs;
%interpolate hydrostatic parameters based on doy
coeffs = avg - amp*cos(doyp);
%repeated fraction expansion of hydrostatic coefficients slant factor
denom = sin(sEL) + coeffs(1)/(sin(sEL) + coeffs(2)/(sin(sEL) + coeffs(3)));
num = 1 + coeffs(1)/(1 + coeffs(2)/(1 + coeffs(3)));
mEh = num / denom;
%repeated fraction expansion of wet coefficients slant factor
denom = sin(sEL) + wet(1)/(sin(sEL) + wet(2)/(sin(sEL) + wet(3)));
num = 1 + wet(1)/(1 + wet(2)/(1 + wet(3)));
mEw = num / denom;
%repeated fraction expansion of height above geoid hydrostatic correction
denom = sin(sEL) + ht(1)/(sin(sEL) + ht(2)/(sin(sEL) + ht(3)));
num = 1 + ht(1)/(1 + ht(2)/(1 + ht(3)));
%note: this should be height above sea level, not height above ellipsoid
%however, it's always within 100m, so it's inside Neill's error budget
dmEh = (1/sin(sEL) - num/denom)*(alt/1000);
mEh = mEh + dmEh;
%calculate tropospheric delays
delTauNa = (znd_hydro*mEh + znd_wet*mEw) / cLight;
