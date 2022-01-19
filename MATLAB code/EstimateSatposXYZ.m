function [Output] = EstimateSatposXYZ(t,PRN_Num,Eph,PRN_Indx,Date,Pn,GPScons)
%% =========================================
% Objective: To estimate the satellite position and satellite clock bias by using the ephemeris data.
% Exemple:[Output] = EstimateSatposXYZ(t,PRN_Num,Eph,PRN_Indx,Date,Pn,GPScons).
% Input:
% t is a considered time (Second of day),
% PRN_Num is a PRN number of the GPS satellites,
% Eph is an ephemeris data,
% PRN_Indx is a considered PRN number,
% Date is a considered date,
% Pn is a pseudorange between antennas of satellite and receiver,
% GPScons is a structure data type of the constant value.
% Output is a structure data type of the satellite position and satellite clock bias.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% =========================================
% Constant
GM = GPScons.GM;	% Earth gravitational constant(m^3/s^2)
We = GPScons.we;	% Earth rotation rate(rad/sec)
c  = GPScons.c;     % Speed of light (m/s)
f1 = GPScons.f1;    % L1 frequency
f2 = GPScons.f2;    % L2 frequency
Y = Date(1);    % Year
M = Date(2);    % Month
D = Date(3);    % Day

% Computed the second of the GPS Week (SOW)
SOW = (day(datetime(Y,M,D),'dayofweek')-1)*(24*60*60);

% Estimated the time to compute the satellite position
Tr = Pn/c;
Ts = mod(SOW + t - Tr,7*86400); % Second of week unit

% Finding the accurate ephemeris
PRN = find(PRN_Indx == PRN_Num);
Toe = Eph(PRN,18);	% Time of Ephemeris(SOW: sec of GPS week)
[~,Row] = min(abs(Ts - Toe));
PRN = PRN(Row);

% === Read Ephemeris data
% Orbit Parameters
a       = Eph(PRN,17).^2;      % Semi-major axis                       (m)
e       = Eph(PRN,15);         % Eccentricity
w0      = Eph(PRN,24);         % Argument of perigee                   (rad)
W0      = Eph(PRN,20);         % Right ascension of ascending node     (rad)
Wdot    = Eph(PRN,25);         % Rate of right ascension               (rad/sec)
i0      = Eph(PRN,22);         % Inclination                           (rad)
idot    = Eph(PRN,26);         % Rate of inclination                   (rad/sec)
M0      = Eph(PRN,13);         % Mean anomaly                          (rad)
delta_n = Eph(PRN,12);         % Mean motion rate                      (rad/sec)

% Corrected coefficients
Cuc     = Eph(PRN,14);         % Argument of perigee (cos)   (rad)
Cus     = Eph(PRN,16);         % Argument of perigee (sine)  (rad)
Crc     = Eph(PRN,23);         % Orbit radius        (cos)   (m)
Crs     = Eph(PRN,11);         % Orbit radius        (sine)  (m)
Cic     = Eph(PRN,19);         % Inclination         (cos)   (rad)
Cis     = Eph(PRN,21);         % Inclination         (sine)  (rad)

% Time
Toe     = Eph(PRN,18);         % Time of Ephemeris(SOW: sec of GPS week)
% GPSW    = Eph(PRN,28);         % GPS Week
% Ttr     = Eph(PRN,34);         % Transmission time of message -604800(SOW:sec of GPS week)
% y       = Eph(PRN,1);          % Year
% m       = Eph(PRN,2);          % Month
% d       = Eph(PRN,3);          % Day of month
% Hour    = Eph(PRN,4);          % Hour
% Min     = Eph(PRN,5);          % Minute
% Sec     = Eph(PRN,6);          % Second

% Clock
af0     = Eph(PRN,7);          % SV Clock Bias(sec)
af1     = Eph(PRN,8);          % SV Clock Drift(sec/sec)
af2     = Eph(PRN,9);          % SV Clock Drift rate(sec/sec^2)
Tgd     = Eph(PRN,32);         % SV Group delay(sec)

% Status
% SV_health   = Eph(PRN,31);     % SV Health
% SV_accuracy = Eph(PRN,30);     % SV Accuracy
% L2_P_flag   = Eph(PRN,29);     % L2 P data flag
% L2_code     = Eph(PRN,27);     % Code on L2 channel
% IODC        = Eph(PRN,33);     % Issue of Data, Clock
% IODE        = Eph(PRN,10);     % Issue of Data, Ephemeris

Tk = Ts - Toe;	% Time elaped since Toe(SOW : sec of GPS week)
M  = M0 + (sqrt(GM/a^3)+delta_n)*Tk;	% Mean anomaly at Tk

% Iterative solution for E
E_old = M;
dE = 1;
while (dE > 10^-12)
    E = E_old - (E_old - e*sin(E_old) - M)/(1 - e*cos(E_old));    % Eccentris anomaly
    dE = abs(E - E_old);
    E_old = E;
end
v = atan2(sqrt(1-e^2)*sin(E),cos(E)-e);	% True anomaly
p = v + w0;	% Argument of perigee + True anomaly or Argument of latitude

% Corrected the orbital errors
u = p + Cuc*cos(2*p) + Cus*sin(2*p);	% Corrected argument of latitude(rad)
r = a*(1-e*cos(E)) + Crc*cos(2*p) + Crs*sin(2*p);	% Corrected radius(meter)
i = i0 + idot*Tk + Cic*cos(2*p) + Cis*sin(2*p);	% Corrected inclination(rad)

W = W0 + (Wdot-We)*Tk - (We*Toe);	% Argument of ascending node

% rotation matrix
R = [cos(W)*cos(u)-sin(W)*sin(u)*cos(i);
    sin(W)*cos(u)+cos(W)*sin(u)*cos(i);
    sin(u)*sin(i)];

Output.XYZ = r.*R;   % Satellite position vector (ECEF with WGS-84)

% === Computed the clock errors
F = -4.442807633e-10; % Constant  (sec/meter^1/2)
Delta_tr = F*e*sqrt(a)*sin(E); % Relativistic correction term
Delta_tsv = af0 + af1*Tk + af2*Tk^2 + Delta_tr - Tgd; % Delta clock
Delta_tsv2 = af0 + af1*Tk + af2*Tk^2 + Delta_tr - (((f1^2)/(f2^2))*Tgd); % Delta clock
dDelta_tsv = af1 + 2*af2*Tk; % Clock drift
Output.t_BiasL1 = Delta_tsv + dDelta_tsv; % The satellite clock bias of the L1 frequency
Output.t_BiasL2 = Delta_tsv2 + dDelta_tsv; % The satellite clock bias of the L2 frequency
end
