function [Out] = PositionA2B(XYZa,XYZb)
% =====================================
% Objective: To compute the elevation and azimuth angles.
% Example: [Out] = PositionA2B(XYZa,XYZb).
% XYZa : Receiver Position in ECEF with WGS-84.
% XYZb : Satellite Positon in ECEF with WGS-84.
% Output is a structure data type of the elevation and azimuth angles.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
% =====================================
LLA = ecef2lla(XYZa(:)');
Lat0  = LLA(1);
Lon0 = LLA(2);
% The rotation matrix R of the ECEF coordinates to the local coordinates
Out.E = [-sind(Lon0) cosd(Lon0) 0];
Out.N = [-sind(Lat0)*cosd(Lon0) -sind(Lat0)*sind(Lon0) cosd(Lat0)];
Out.U = [cosd(Lat0)*cosd(Lon0) cosd(Lat0)*sind(Lon0) sind(Lat0)];

r = [XYZb(1)-XYZa(1) XYZb(2)-XYZa(2) XYZb(3)-XYZa(3)];
Out.P = r./norm(r); % Unit vector

% The ENU => Elevation & Azimuth 
Out.Ele = asin(Out.P*Out.U')*180/pi; % Elevation angle(deg)
Out.Azi = atan2(Out.P*Out.E',Out.P*Out.N')*180/pi; % Azimuth angle(deg)
end
% referenceEllipsoid
