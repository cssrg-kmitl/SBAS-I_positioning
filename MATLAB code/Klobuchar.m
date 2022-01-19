function [Tion] = Klobuchar(Alpha,Beta,Ele,Azi,Lat,Lon,t_gps)
%% =======================================================
% Objective: To estimate the ionospheric delay based on the Klobuchar model.
% Example: [Tion] = Klobuchar(Alpha,Beta,Ele,Azi,Lat,Lon,t_gps).
% Alpha = ION ALPHA
% Beta = ION BETA
% Ele = Elevation angle in degree
% Azi = Azimuth angle in degree
% Lat = latitude in degree
% Long = Longitude in degree
% SOD = Second of day
% Tion = Output
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (February 2019).
%% =======================================================
% === Ionospheric model (20.3.3.5.2.5 of IS-GPS-200D)
El = Ele./180; % Elevation angle (semi-circles)
Az = Azi./180; % Azimuth angle (semi-circles)
Shi = 0.0137./(El + 0.11) - 0.022;  % Earth's central angle
Lam_u = Lon/180;    % User geodetic longitude(semi-circles) WGS-84
Phi_u = Lat/180;      % User geodetic latitude(semi-circles) WGS-84
Phi_i = Phi_u + Shi.*cos(Az.*pi);   % Geodetic latitude (semi-circles)
Phi_i = (Phi_i > 0.416).*0.416 + (Phi_i < -0.416).*(-0.416) +...
    (abs(Phi_i) <= 0.416).*Phi_i;
Lam_i = Lam_u + Shi.*(sin(Az.*pi)./cos(Phi_i.*pi)); % Geodetic longitude (semi-circles)
Phi_m = Phi_i + 0.064.*cos((Lam_i - 1.617).*pi);  % Geomagnetic latitude(assumed height 350 km) (semi-circles)
t_m = (Lam_i.*4.32e4) + t_gps;  % Local time(sec)
t_m = (t_m >= 86400).*(t_m - 86400) + (t_m < 0).*(t_m + 86400) +...
    (t_m < 86400).*(t_m >= 0).*t_m;
F = 1 + 16.*(0.53 - El).^3;    % Obliquity factor
PER = zeros(1,length(Phi_m));
AMP = zeros(1,length(Phi_m));
for i = 1:length(Phi_m)
    PER(i) = Beta(1) + Beta(2)*Phi_m(i) + Beta(3)*Phi_m(i)^2 + Beta(4)*Phi_m(i)^3;
    AMP(i) = Alpha(1) + Alpha(2)*Phi_m(i) + Alpha(3)*Phi_m(i)^2 + Alpha(4)*Phi_m(i)^3;
end
PER = (PER < 72000).*72000 + (PER >= 72000).*PER;
AMP = (AMP >= 0).*AMP;
Xm = 2*pi*(t_m - 50400)./PER;
Tion = (abs(Xm) <= 1.57).*(F.*((5e-9) + AMP.*(1 - (Xm.^2)./2 + (Xm.^4)./24))) +...
    (abs(Xm) > 1.57).*(F.*5e-9);
