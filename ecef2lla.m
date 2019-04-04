function [lat, lon, alt] = ecef2lla(pVec)
% ecef2lla : Convert from a position vector in the Earth-centered,
%   Earth-fixed (ECEF) reference frame to latitude, longitude, and altitude
%   (geodetic with respect to the WGS-84). 
%
%
% INPUTS
%
% pVec ---- 3-by-1 position coordinate vector in the ECEF reference frame,
%           in meters
%
% OUTPUTS
%
% lat ---- latitude in radians
%
% lon ---- longitude in radians
%
% alt ---- altitude (height) in meters above the ellipsoid
%
%
%=======================================================================%

x = pVec(1);
y = pVec(2);
z = pVec(3);

lon = atan2(y, x);
p = sqrt(x.^2 + y.^2);
r = sqrt(x.^2 + y.^2 + z.^2);

a = 6378137; %semi major axis (m)
f = 1/298.257223563; %flattening earth factor
w = 7292115e-11; %mean angular velocity of Earth (rad/sec)
e_sq = 2*f - f^2;

lat = zeros(length(z),1);

for i = 1:length(z)

    lat(i) = asin(z(i)/r(i)); %guess
    latold = inf;

    while abs(lat(i) - latold) > 10^(-12)
        Rn = a/sqrt(1 - e_sq*sin(lat(i))^2);
        latold = lat(i);
        lat(i) = atan2(z(i) + Rn*e_sq*sin(lat(i)), p(i));
    end
end

alt = p./cos(lat) - Rn;

