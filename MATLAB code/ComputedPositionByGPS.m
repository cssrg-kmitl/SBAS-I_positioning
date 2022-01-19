function [Pos] = ComputedPositionByGPS(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,GPScons,K)
%% ================================================
% Objective: To estimate the receiver position by using the GPS information.
% Example: [Pos] = ComputedPositionByGPS(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,GPScons,K).
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (February 2019).
%% ================================================
% === Constant parameters
c  = GPScons.c;             % light speed = 299792458 m/s
Ele_Cut = GPScons.ele_cut;  % Elevation cut-off
Alpha = Nav.Ion.GPSA;       % ION ALPHA
Beta = Nav.Ion.GPSB;        % ION BETA
LLAu = ecef2lla(XYZu(:)');  % Pre-positioning estimation

% === Estimated the user position
XYZs = nan(3,length(PRN));      % Satellite Position
t_Bias = nan(1,length(PRN));	% Satellite Clock Bias
TC = nan(1,length(PRN));        % Tropospheric delays
IC = nan(1,length(PRN));        % Ionospheric delays
t_gps = nan(1,length(PRN));     % GPS time
Ele = nan(1,length(PRN));       % Elevation engle
Azi = nan(1,length(PRN));       % Azimuth engle

for n = 1:length(PRN)
    % Satellite positions (XYZ) and satellite clock bias (sec)
    [sat] = EstimateSatposXYZ(SOD(n),PRN(n),Nav.Eph.G',Nav.PRN.G,Obs.Date.St,Pn(n),GPScons);
    XYZs(1:3,n) = sat.XYZ;
    t_Bias(n) = sat.t_BiasL1;
    [Out] = PositionA2B(XYZu,XYZs(1:3,n));
    Ele(n) = Out.Ele;
    Azi(n) = Out.Azi;
    t_gps(n) = SOD(n) - Pn(n)/c;
end

if K>2
    % Filter with elevation
    Indx_Cut = Ele<Ele_Cut;
    XYZs(:,Indx_Cut) = [];  % Satellite position
    Pn(Indx_Cut) = [];      % Pseudorange
    t_Bias(Indx_Cut) = [];  % Sat clock bias
    t_gps(Indx_Cut) = [];   % Sat clock bias
    TC(Indx_Cut) = [];      % Tropospheric delays
    IC(Indx_Cut) = [];      % Ionospheric delays
    Ele(Indx_Cut) = [];     % Elevation angle
    Azi(Indx_Cut) = [];     % Elevation angle
    PRN(Indx_Cut) = [];     % PRN numbar
end

if (length(PRN) > 3) % Check number of satellite
    if K>2
        % Estimated the ionospheric corrections (20.3.3.5.2.5 of IS-GPS-200D)
        [Tion] = Klobuchar(Alpha,Beta,Ele,Azi,LLAu(1),LLAu(2),t_gps);
        IC = -c*Tion;
        % Estimated the Tropospheric corrections (A.4.2.4 of RTCA DO-229D)
        [d_hyd,d_wet] = EstimateTropDalay(LLAu(1),Height,DOY);
        TC = -(d_hyd + d_wet)*(1.001./sqrt(0.002001 + sind(Ele).^2)); % Tropospheric model(A.4.2.4 of RTCA DO-229D)
    else
        IC(:) = 0;
        TC(:) = 0;
    end

    % Pseudorange
    R = Pn + c.*(t_Bias) + TC + IC;
    
    % Sagnac effect correction
    zXYZs = nan(3,length(XYZs));
    [~,n] = size(XYZs);
    for i = 1:n
        Time_Rota = R(i)/c;
        zXYZs(:,i) = Rotation_Z(XYZs(:,i),Time_Rota);
    end
    
    cIter = 0;      % Interative count
    nIter = 100;    % Iterative number
    dX = inf;       % Initial dX
    Bias = 0;       % Initial bias
    while (norm(dX) > 10^-3) && (cIter < nIter) % Iterative method
        cIter = cIter + 1;
        
        % Computed the line-of-sight vectors and ranges from satellite to receiver
        pj = zXYZs - XYZu(:)*ones(1,size(zXYZs,2));
        Norm_pj = sqrt(sum(pj.^2,1));
        abc = -pj./(ones(3,1)*Norm_pj); % line of sight unit vectors
        
        % Computed the a-priori range
        P_hat = Norm_pj' + Bias;
        P = R' - P_hat;                  % matrix P
        
        A = [abc',ones(length(abc),1)]; % matrix A
        
        % MMSE solution
        dX = (A'*A)\(A'*P);
        
        % Update XYZu and Bias
        XYZu = XYZu(:) + dX(1:3);
        Bias = Bias + dX(4);
    end
    
    LLAu = ecef2lla(XYZu'); % User position in the geodetic coordinate
    Pos.Sat_Num = sum(PRN>0);
    Pos.PRN = PRN;
    Pos.llh = LLAu;
    Pos.xyz = XYZu;
    Pos.bias = Bias;
else
    Pos = [];
end
