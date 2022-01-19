function [Pos] = ComputedPositionByGPSwithSBAS_FL(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,SBAS,GPScons,K)
%% ================================================
% Objective: To estimate the receiver position by using the GPS information with the SBAS corrections of the fast and long-term corrections.
% Example: [Pos] = ComputedPositionByGPSwithSBAS_FL(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,SBAS,GPScons,K).
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (February 2019).
%% ================================================
% === Constant parameters
c  = GPScons.c;             % light speed (299792458 m/s)
Ele_Cut = GPScons.ele_cut;  % Elevation cut-off
Alpha = Nav.Ion.GPSA;       % ION ALPHA
Beta = Nav.Ion.GPSB;        % ION BETA
LLAu = ecef2lla(XYZu(:)');  % Pre-positioning estimation

% === SBAS parameters
Time_pr = SBAS.Time_pr; % Time of applicability of a previous fast corrections
Time_ap = SBAS.Time_ap; % Time of applicability of the most recent fast corrections
oPRC = SBAS.oPRC;
PRC = SBAS.PRC;
Sig2_UDRE = SBAS.Sig2_UDRE;
Indx_NAN = isnan(Sig2_UDRE(PRN));
PRN(Indx_NAN) = [];
Pn(Indx_NAN) = [];
SOD(Indx_NAN) = [];
t_lat = SBAS.t_lat;
ai2 = SBAS.ai2;
Ifc1 = SBAS.Ifc1;
Ifc2 = SBAS.Ifc2;
Brre = SBAS.Brre;
Cltc_1sb = SBAS.Cltc_1sb;
Cltc_v1 = SBAS.Cltc_v1;
Iltc_v1 = SBAS.Iltc_v1;
Cltc_v0 = SBAS.Cltc_v0;
Iltc_v0 = SBAS.Iltc_v0;
Cgeo_1sb = SBAS.Cgeo_1sb;
Cgeo_v = SBAS.Cgeo_v;
Igeo = SBAS.Igeo;
Cer = SBAS.Cer;
Ciono_step = SBAS.Ciono_step;
Iiono = SBAS.Iiono;
Ciono_ramp = SBAS.Ciono_ramp;
RSSu = SBAS.RSSu;
RSSi = SBAS.RSSi;
C_cov = SBAS.C_cov;
Vcod = SBAS.Vcod;
d_XYZ = SBAS.d_XYZ;
d_XYZr = SBAS.d_XYZr;
d_af0 = SBAS.d_af0;
d_af1 = SBAS.d_af1;
t_o = SBAS.t_o;
for o = PRN(:)'
    % Eliminated the empty values
    if sum(isnan(d_af0(o)))>0
        PRC(o)= 0;
        d_XYZ{1,o} = zeros(3,1);
        d_XYZr{1,o} = zeros(3,1);
        d_af0(o) = 0;
        d_af1(o) = 0;
        t_o(o) = 0;
    end
end
Cov = SBAS.Cov;
SF = SBAS.SF;

% === Estimated the user position
XYZs = nan(3,length(PRN));      % Satellite Position
Long_XYZ = nan(3,length(PRN));	% Long-term corrections from SBAS
t_Bias = nan(1,length(PRN));    % Satellite Clock Bias
d_t = nan(1,length(PRN));       % Satellite Clock Bias corrections from SBAS
T_Bias = nan(1,length(PRN));	% Clock time errors
TC = nan(1,length(PRN));        % Tropospheric delays
IC = nan(1,length(PRN));        % Ionospheric delays
t_gps = nan(1,length(PRN));     % GPS time
Ele = nan(1,length(PRN));       % Elevation engle
Azi = nan(1,length(PRN));       % Azimuth engle

for n = 1:length(PRN)
    % Satellite positions (XYZ) and satellite clock bias (sec)
    [Sat] = EstimateSatposXYZ(SOD(n),PRN(n),Nav.Eph.G',Nav.PRN.G,Obs.Date.St,Pn(n),GPScons);
    Long_XYZ(1:3,n) = d_XYZ{PRN(n)} + d_XYZr{PRN(n)}.*ones(3,1)*(SOD(n) - t_o(PRN(n)));
    d_t(n) = d_af0(PRN(n)) + d_af1(PRN(n))*(SOD(n) - t_o(PRN(n)));
    XYZs(1:3,n) = Sat.XYZ;
    t_Bias(n) = Sat.t_BiasL1;
    T_Bias(n) = t_Bias(n) + d_t(n);
    [Out] = PositionA2B(XYZu,XYZs(1:3,n));
    Ele(n) = Out.Ele;
    Azi(n) = Out.Azi;
    t_gps(n) = SOD(n) - Pn(n)/c;
end

if K > 2
    % Estimated the ionospheric delays (20.3.3.5.2.5 of IS-GPS-200D)
    [Tion] = Klobuchar(Alpha,Beta,Ele,Azi,LLAu(1),LLAu(2),t_gps); 
    IC = -c*Tion;
    % Estimated the Tropospheric delays (A.4.2.4 of RTCA DO-229D)
    [d_hyd,d_wet] = EstimateTropDalay(LLAu(1),Height,DOY);
    TC = -(d_hyd + d_wet)*(1.001./sqrt(0.002001 + sind(Ele).^2));
else
    IC(:) = 0;
    TC(:) = 0;
end

RRC = (PRC(PRN) - oPRC(PRN))./(Time_ap(PRN) - Time_pr(PRN)); % Range-rate corrections (A.4.4.3 of RTCA DO-229D)
dT = (SOD(:)' - Time_ap(PRN)); % Differentail time (A.4.4.3 of RTCA DO-229D)

% Corrected pseudoranges
R = Pn + c.*(T_Bias) + TC + IC + PRC(PRN) + RRC.*dT;

% Sagnac effect corrections
zXYZs = nan(3,length(XYZs));
[~,n] = size(XYZs);
for i = 1:n
    Time_Rota = R(i)/c;
    zXYZs(:,i) = Rotation_Z(XYZs(:,i)+Long_XYZ(:,i),Time_Rota);
end

Weight = 1./Sig2_UDRE(PRN); % Weighting matrix

if K > 2
    % Filter with elevations
    Indx_Cut = Ele<Ele_Cut;
    zXYZs(:,Indx_Cut) = [];
    R(Indx_Cut) = [];
    Weight(Indx_Cut) = [];
end

if (length(R) > 3) % Check number of satellite
    cIter = 0;      % Interative counter
    nIter = 100;    % Iterative number
    dX = inf;       % Initial dX
    Bias = 0;       % Initial bias
    W = diag(Weight);
    while (norm(dX) > 10^-3) && (cIter < nIter) % Iterative method
        cIter = cIter + 1;
        
        % Computed line-of-sight vectors and ranges from satellite to receiver
        pj = zXYZs - XYZu(:)*ones(1,size(zXYZs,2));
        Norm_pj = sqrt(sum(pj.^2,1));
        abc = -pj./(ones(3,1)*Norm_pj); % Lin-of-sight unit vectors
        
        % Computed the a-priori range
        P_hat = Norm_pj' + Bias;
        P = R' - P_hat;                  % matrix P
        
        A = [abc',ones(length(abc),1)]; % matrix A
        
        % MMSE solution
        dX = (A'*W*A)\(A'*W*P);
                
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
