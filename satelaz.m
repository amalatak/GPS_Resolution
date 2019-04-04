function [elRad, azRad] = satelaz(rSvEcef,rRxEcef)
% satelaz : Compute satellite elevation and azimuth angles in radians with
% respect to receiver location.
%
%
% INPUTS
%
% rSvEcef ---- 3-by-1 satellite location in ECEF coordinates, in meters.
%
% rRxEcef ---- 3-by-1 receiver location in ECEF coordinates, in meters.
%
%
% OUTPUTS
%
% elRad ------ Satellite elevation angle (the angle between the WGS-84 local
% ENU tangent plane and the receiver-to-satellite vector),
% in radians.
%
% azRad ------ Satellite azimuth angle (the angle in the ENU tangent plane
% between North and the receiver-to-satellite vector projection,
% measured positive clockwise), in radians.
%
%+------------------------------------------------------------------------------+
% References:
% Moriba Jah, Ph.D Lecture Notes, ASE 372N Fall 2018
%
% Author: Andrew Malatak
%+==============================================================================+

% receiver lat and long in radians and meters:

[latr, lonr, ~] = ecef2lla(rRxEcef);


% enu up vector

e_enu_up = [cos(latr)*cos(lonr); cos(latr)*sin(lonr); sin(latr)];


% ecef to enu rotation matrix 

Recef2enu = ecef2enu(latr, lonr);

% line of sight unit vector
if size(rSvEcef) == size(rRxEcef)
	e_LOS = (rSvEcef - rRxEcef)./norm(rSvEcef - rRxEcef);
else
	e_LOS = (rSvEcef - rRxEcef')./norm(rSvEcef' - rRxEcef');
end

% find enu unit vector and then az and el
e_enu = Recef2enu*e_LOS;

azRad = atan2(e_enu(1), e_enu(2));
elRad = asin(dot(e_enu_up, e_LOS));

