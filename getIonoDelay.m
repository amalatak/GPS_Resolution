function [delTauG] = getIonoDelay(ionodata,fc,rRx,rSv,tGPS,model)
% getIonoDelay : Return a model-based estimate of the ionospheric delay
%                experienced by a transionospheric GNSS signal as it
%                propagates from a GNSS SV to the antenna of a terrestrial
%                GNSS receiver.
%
% INPUTS
%
% ionodata ------- Structure containing a parameterization of the
%                  ionosphere that is valid at time tGPS.  The structure is
%                  defined differently depending on what ionospheric model
%                  is selected:
%
%                  broadcast --- For the broadcast (Klobuchar) model, ionodata
%                  is a structure containing the following fields: 
% 
%                     alpha0 ... alpha3 -- power series expansion coefficients
%                     for amplitude of ionospheric TEC
%
%                     beta0 .. beta3 -- power series expansion coefficients
%                     for period of ionospheric plasma density cycle
%
%
%                  Other models TBD ...
%
% fc ------------- Carrier frequency of the GNSS signal, in Hz.
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
%                  broadcast --- The broadcast (Klobuchar) model.  
%
%                  Other models TBD ...
%
% OUTPUTS
%
% delTauG -------- Modeled scalar excess group ionospheric delay experienced
%                  by the transionospheric GNSS signal, in seconds.
% 
%+------------------------------------------------------------------------------+
% References: For the broadcast (Klobuchar) model, see IS-GPS-200F
% pp. 128-130.
%
%+==============================================================================+
% Copyright © 2006 Isaac Miller, all rights reserved.  No portion of this code
% may be reproduced, reused, or distributed in any form without prior written
% permission of the author.
  
if(strcmp(model,'broadcast')~=1)
  error('Unrecognized model.');
end
  
% extract ionospheric parameters
alpha0 = ionodata.alpha0;
alpha1 = ionodata.alpha1;
alpha2 = ionodata.alpha2;
alpha3 = ionodata.alpha3;
beta0 = ionodata.beta0;
beta1 = ionodata.beta1;
beta2 = ionodata.beta2;
beta3 = ionodata.beta3;
% receiver geodetic latitude and longitude, in radians
[lat,lon,alt] = ecef2lla(rRx);
% satellite elevation and azimuth
[sEL, sAZ] = satelaz(rSv,rRx);
% earth's central angle between receiver and earth projection of ionospheric
% intersection point, in radians
psi = 0.0137*pi/(sEL/pi + 0.11) - 0.022*pi;
% slant factor to account for non-straight path thru ionosphere
Fs = 1 + 16*(0.53 - sEL/pi)^3;
% ionospheric intersection point geodetic latitude (rad.)
phi_i = lat + psi*cos(sAZ);
if (abs(phi_i) > pi*0.416)
    phi_i = pi*0.416*sign(phi_i);
end
% ionospheric intersection point geodetic longitude (rad.)
lamda_i = lon + psi*sin(sAZ)/cos(phi_i);
% geomagnetic latitude (rad.)
phi_m = phi_i + 0.064*pi*cos(lamda_i - 1.617*pi);
% geomagnetic latitude (semicircles) for expansion coefficients
phi_m_ss = phi_m / pi;
% local time of day in seconds (local offset + time at Greenwich)
tlocal = 43200*lamda_i/pi + mod(tGPS.seconds,86400);
if (tlocal < 0)
    tlocal = tlocal + 86400;
end
if (tlocal >= 86400)
    tlocal = tlocal - 86400;
end
% ionospheric period (sec)
T = beta0 + beta1*phi_m_ss + beta2*phi_m_ss^2 + beta3*phi_m_ss^3;
if (T < 72000)
    T = 72000;
end
% local time angle (rad.)
theta = 2*pi*(tlocal - 50400)/T;
% amplitude of ionospheric fluctuation
C1 = alpha0 + alpha1*phi_m_ss + alpha2*phi_m_ss^2 + alpha3*phi_m_ss^3;
if (C1 < 0)
    C1 = 0;
end
% final ionospheric delay (sec)
if (abs(theta) < 1.57)
    delTauG = Fs * (5e-9 + C1*(1 - theta^2/2 + theta^4/24));
else
    delTauG = Fs * 5e-9;
end

% The above calculations are for the GPS L1 frequency.  Multiply delTauG by
% (fL1^2)/(fc^2) to map to the carrier frequency fc.
fL1 = 1575.42e6;
delTauG = delTauG*((fL1/fc)^2);