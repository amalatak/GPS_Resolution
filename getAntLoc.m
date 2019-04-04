function [rEcef,vEcef,epochYear] = getAntLoc(antId)
% getAntLoc : Get the location of a GNSS antenna from a list of known
%             antennas.
%
% INPUTS
%
% antId ------ Antenna identifier.  Must be among those in the list below.
%
%
% OUTPUTS
%
% rEcef ------ 3-by-1 cartesian position vector corresponding to the L1
%              antenna phase center in the ITRF08 reference frame at
%              epochYear.
%
% vEcef ------ 3-by-1 cartesian velocity vector corresponding to the L1
%              antenna phase center relative to the ITRF08 reference frame
%              and expressed in the ITRF08 reference frame, in meters/yr.
%              Propagate rEcef to the current epoch by
%
%              rEcefCurrent = rEcef + vEcef*(currentYear - epochYear)
%
% epochYear -- 4-digit year and decimal fraction of year of epoch at which
%              rEcef is defined (e.g., 2005.0). 
% 
% 
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys
%+==============================================================================+
    
antId = upper(antId);

switch antId
  case 'TXAU'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=TXAU
    rEcef = [-743774.390; -5460644.512; 3200347.710];
    vEcef = [-0.0120;  -0.0003; -0.0031];
    epochYear = 2005.0;
  case 'SAM2'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=SAM2
    rEcef = [-751970.424; -5463638.752; 3193399.042];
    vEcef = [-0.0123; -0.0003; -0.0033];
    epochYear = 2005.0;
  case 'TXTA'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates: 
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=TXTA
    rEcef = [-712254.169; -5450502.920; 3224462.604];
    vEcef = [-0.0114; -0.0024; -0.0027];
    epochYear = 2005.0;
  case 'WRW0'
    % WRW primary antenna. Type: TRM57971.00. Source of coordinates: L1-only CDGPS
    % solution on the WRW0-TXAU baseline using the above TXAU coordinates.
    % Velocity copied from TXAU velocity.
    rEcef = [-741990.378;  -5462227.718;  3198019.532];
    vEcef = [-0.0120;  -0.0003; -0.0031];
    epochYear = 2005.0;
  case 'WRW1'
    error('Not yet defined');
  case 'AUT0'  
    % CSR primary antenna. Type: TRM57971.00.  Source of coordinates: Integrated
    % MGEX network solution.
    rEcef = [-741133.4180; -5456322.9775; 3208395.3863];
    vEcef = [-0.0120;  -0.0003; -0.0031];
    epochYear = 2005.0;
  case 'ARL0'
    % ARL primary antenna. Type: TRM57971.00.  Source of coordinates: L1-only
    % CDGPS solution on the WRW0-ARL0 baseline using the above WRW0
    % coordinates.
    rEcef = [-740314.629; -5457071.504; 3207240.3226];
    vEcef = [-0.0120;  -0.0003; -0.0031];
    epochYear = 2005.0;
  case 'RHO0'
    % Rhodes north antenna.  Rough CASES solution. 
    rEcef = [1101936.511; -4583477.965; 4282260.433];
    vEcef = [0;0;0];
    epochYear = 2005.0;
  case 'RHO1'
    % Rhodes south antenna. Rough CASES solution. 
    rEcef = [1101968.234; -4583484.254; 4282240.143];
    vEcef = [0;0;0];
    epochYear = 2005.0;
  case 'PHI0'
    % Philips northeast antenna. Rough CASES solution. 
    rEcef = [1101898.959; -4583390.505; 4282331.985];
    vEcef = [0;0;0];
    epochYear = 2005.0;
  case 'FICTION_ISLAND'
    rEcef = [6378137; 0; 0];
    vEcef = [0;0;0];
    epochYear = 2005.0;
  case 'ARP7'
    % NGS CORS antenna. Type: TRM41249USCG.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=ARP7
    rEcef = [-693606.128; -5601312.014; 2960669.097];
    vEcef = [-0.0137; 0.0004; -0.0037];
    epochYear = 2005.0;
  case 'ARP8'
    % NGS CORS antenna. Type: TRM41249USCG.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=ARP8
    rEcef = [-693632.879; -5601307.978; 2960670.524];
    vEcef = [-0.0127; 0.0013; -0.0028];
    epochYear = 2005.0;
  case 'TXPO'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=TXPO
    rEcef = [-694669.291; -5601118.520; 2960775.720];
    vEcef = [-0.0116; -0.0003; -0.0032];
    epochYear = 2005.0;
  case 'TXCC'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=TXPO
    rEcef = [-731659.034; -5601556.959; 2951108.120];
    vEcef = [-0.0122; 0.0028; -0.0030];
    epochYear = 2005.0;
  case 'TXRP'
    % NGS CORS antenna. Type: TRM57971.00.  Source of coordinates:
    % http://www.ngs.noaa.gov/cgi-cors/corsage.prl?site=TXPO
    rEcef = [-691206.427; -5589900.949; 2982564.951];
    vEcef = [-0.0130; -0.0001; -0.0028];
    epochYear = 2005.0;
  case 'EH20'
    % From flyday survey on the East Hatch at the 20 yard line.  See data in
    % /vtrak/data2/flydayJul25/surveyEH20.  Source of coordinates: L1-only
    % CDGPS solution on the WRW0-EH20 baseline using the above WRW0
    % coordinates.  Velocity copied from TXAU velocity.
    rEcef = [-741710.846; -5462470.21628; 3197610.8423];
    vEcef = [-0.0120;  -0.0003; -0.0031];
    epochYear = 2005.0;
  otherwise
    error('Unrecognized antenna identifier');
end
