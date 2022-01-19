function [Pos,SD] = ComputedPositionByGPSwithSBAS_FLI(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,SBAS,GPScons,K)
%% ================================================
% Objective: To estimate the receiver position by using the GPS information with the SBAS corrections of the fast and
% long-term corrections + Ionospheric grid point.
% Example: [Pos,SD] = ComputedPositionByGPSwithSBAS_FLI(DOY,SOD,PRN,Obs,Nav,Pn,XYZu,Height,SBAS,GPScons,K).
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (February 2019).
%% ================================================
% === Constant parameters
c  = GPScons.c;             % light speed (299792458 m/s)
Re = GPScons.Rex;           % The radius of the Earth (6378136.3 m)
h = GPScons.h;              % Ionospheric height for mapping function (350 km)
Ele_Cut = GPScons.ele_cut;  % Elevation cut-off
LLAu = ecef2lla(XYZu(:)');  % Pre-positioning estimation

% === SBAS parameters ===
Time_pr = SBAS.Time_pr;
Time_ap = SBAS.Time_ap;
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
to = SBAS.t_o;
for o = PRN(:)'
    if sum(isnan(d_af0(o)))>0
        PRC(o)= 0;
        d_XYZ{1,o} = zeros(3,1);
        d_XYZr{1,o} = zeros(3,1);
        d_af0(o) = 0;
        d_af1(o) = 0;
        to(o) = 0;
    end
end
LAT = SBAS.Lat;
LON = SBAS.Lon;
IVDE = SBAS.IVDE;
Sig2_GIVE = SBAS.Sig2_GIVE;
Cov = SBAS.Cov;
SF = SBAS.SF;

% === Calculated the positioning correction ===
% prepare matrix
XYZs = nan(3,length(PRN));      % Satellite Position
Long_XYZ = nan(3,length(PRN));	% Satellite Position Long-term
t_Bias = nan(1,length(PRN));	% Satellite Clock Bias
d_t = nan(1,length(PRN));       % Satellite Clock Bias from SBAS
T_Bias = nan(1,length(PRN));	% Clock time errors
TC = nan(1,length(PRN));        % Tropospheric delays
IC = nan(1,length(PRN));        % Ionospheric delay
Sig2_UIRE = nan(1,length(PRN));	% Standard deviation of Ionospheric delay
t_gps = nan(1,length(PRN));     % GPS time
Ele = nan(1,length(PRN));       % Elevation engle
Azi = nan(1,length(PRN));       % Azimuth engle

for n = 1:length(PRN)
    % Satpos XYZ and satellite clock bias (sec)
    [Sat] = EstimateSatposXYZ(SOD(n),PRN(n),Nav.Eph.G',Nav.PRN.G,Obs.Date.St,Pn(n),GPScons);
    Long_XYZ(1:3,n) = d_XYZ{PRN(n)} + d_XYZr{PRN(n)}.*ones(3,1)*(SOD(n) - to(PRN(n)));
    d_t(n) = d_af0(PRN(n)) + d_af1(PRN(n))*(SOD(n) - to(PRN(n)));
    XYZs(1:3,n) = Sat.XYZ;
    t_Bias(n) = Sat.t_BiasL1;
    T_Bias(n) = t_Bias(n) + d_t(n);
    [Out] = PositionA2B(XYZu,XYZs(1:3,n));
    Ele(n) = Out.Ele;
    Azi(n) = Out.Azi;
    t_gps(n) = SOD(n) - Pn(n)/c;
end

if K > 2
    Fpp = (1 - (Re.*cosd(Ele)./(Re + h)).^2).^-0.5;   % Obliquity factor (A.4.4.10.4 of RTCA DO-229D)
    % Ionospheric delays from IGP based on the interpolated method
    Var_Lat = [];
    Var_Lon = [];
    Var_IVDE = [];
    Var_Sig2_GIVE = [];
    % Find Band
    LLAu(mod(LLAu,1)==0) = LLAu(mod(LLAu,1)==0)+0.000001;
    if LLAu(2)>-180.1 && LLAu(2)<=-140 && abs(LLAu(1))<55
        Band = [9 1 2];
    elseif LLAu(2)>-140 && LLAu(2)<=-100 && abs(LLAu(1))<55
        Band = [1 2 3];
    elseif LLAu(2)>-100 && LLAu(2)<=-60 && abs(LLAu(1))<55
        Band = [2 3 4];
    elseif LLAu(2)>-60 && LLAu(2)<=-20 && abs(LLAu(1))<55
        Band = [3 4 5];
    elseif LLAu(2)>-20 && LLAu(2)<=20 && abs(LLAu(1))<55
        Band = [4 5 6];
    elseif LLAu(2)>20 && LLAu(2)<=60 && abs(LLAu(1))<55
        Band = [5 6 7];
    elseif LLAu(2)>60 && LLAu(2)<=100 && abs(LLAu(1))<55
        Band = [6 7 8];
    elseif LLAu(2)>100 && LLAu(2)<=140 && abs(LLAu(1))<55
        Band = [7 8 9];
    elseif LLAu(2)>140 && LLAu(2)<=180 && abs(LLAu(1))<55
        Band = [8 9 1];
    end
    for i = Band
        Indx = find(~isnan(LAT{1,i}));
        if ~isempty(Indx)
            Var_Lat = [Var_Lat LAT{1,i}(Indx)];
            Var_Lon = [Var_Lon LON{1,i}(Indx)];
            Var_IVDE = [Var_IVDE IVDE{1,i}(Indx)];
            Var_Sig2_GIVE = [Var_Sig2_GIVE Sig2_GIVE{1,i}(Indx)];
        end
    end
    Indx_Lat = abs(Var_Lat-LLAu(1)) < 5 & abs(Var_Lat-LLAu(1)) >= 0;
    Indx_Lon = abs(Var_Lon-LLAu(2)) < 5 & abs(Var_Lon-LLAu(2)) >= 0;
    Indx = Indx_Lat&Indx_Lon;
    LAT = Var_Lat(Indx);
    LON = Var_Lon(Indx);
    Buff_IVDE = Var_IVDE(Indx);
    Buff_uive = Var_Sig2_GIVE(Indx);
    % Interpolatiom algorithm
    if length(Buff_IVDE)==3 % Three-points weight
        Num_LAT = histcounts(LAT,2);
        Num_LON = histcounts(LON,2);
        cLAT = unique(LAT);
        cLON = unique(LON);
        [~,LAT_1] = max(Num_LAT);
        [~,LAT_2] = min(Num_LAT);
        [~,LON_1] = max(Num_LON);
        [~,LON_2] = min(Num_LON);
        Ypp = (LLAu(1)-cLAT(LAT_1))./(cLAT(LAT_2)-cLAT(LAT_1));
        Xpp = (LLAu(2)-cLON(LON_1))./(cLON(LON_2)-cLON(LON_1));
        W_ion = [Ypp 1-Xpp-Ypp Xpp];
        IC_vpp = W_ion*Buff_IVDE';
        Sig2_uive = W_ion*Buff_uive';
    elseif length(Buff_IVDE)==4 % Four-points weight
        Ypp = (LLAu(1)-LAT(1))./(LAT(2)-LAT(1));
        Xpp = (LLAu(2)-LON(1))./(LON(3)-LON(1));
        W_ion = [(1-Xpp)*(1-Ypp) (1-Xpp)*Ypp (1-Ypp)*Xpp  Ypp*Xpp];
        IC_vpp = W_ion*Buff_IVDE';
        Sig2_uive = W_ion*Buff_uive';
    else % Default
        IC_vpp = Buff_IVDE;
        Sig2_uive = Buff_uive;
    end
    
    if isnan(IC_vpp) % Nan case
        IC_vpp = 0;
        Sig2_uive = 0;
    elseif isempty(IC_vpp) % Empty case
        IC_vpp = 0;
        Sig2_uive = 0;
    end
    
    IC = -Fpp.*IC_vpp;
    Sig2_UIRE = (Fpp.^2).*Sig2_uive;
    
    % Estimated the Tropospheric delays (A.4.2.4 of RTCA DO-229D)
    [d_hyd,d_wet] = EstimateTropDalay(LLAu(1),Height,DOY);
    TC = -(d_hyd + d_wet)*(1.001./sqrt(0.002001 + sind(Ele).^2)); % Tropospheric model(A.4.2.4 of RTCA DO-229D)
else
    IC(:) = 0;
    TC(:) = 0;
    Sig2_UIRE(:) = 0;
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

% Toltal sigma^2
if K > 2
    XYZ_R = zXYZs - XYZu(:)*ones(1,size(XYZs,2)); % Pseudorange vector between antennas of satellite and receiver
    Norm_XYZ = sqrt(sum(XYZ_R.^2,1));
    Unit_XYZ = XYZ_R./(ones(3,1)*Norm_XYZ); % line of sight unit vectors
    I = [Unit_XYZ; ones(1,size(Unit_XYZ,2))];
    if RSSu == 1
        Ec = C_cov.*SF(PRN);
        d_UDRE = nan(1,length(PRN));
        Efc = ai2(PRN).*(dT - t_lat(PRN).^2)./2;
        Errc = zeros(1,length(PRN));
        Eltc = zeros(1,length(PRN));
        Eer = zeros(1,length(PRN));
        for o = 1:length(PRN)
            d_UDRE(o) = sqrt(I(:,o)'*Cov{1,PRN(o)}*I(:,o)) + Ec(o);
        end
        d_UDRE(isnan(d_UDRE)) = 1;
        Sig2_flt = (sqrt(Sig2_UDRE(PRN)).*d_UDRE).^2 + Efc.^2 + Errc.^2 + Eltc.^2 + Eer.^2; % (A.4.5.1 of RTCA DO-229D)
    elseif RSSu == 0
        Ec = C_cov.*SF(PRN);
        d_UDRE = nan(1,length(PRN));
        Efc = ai2(PRN).*(dT - t_lat(PRN).^2)./2;
        Errc = zeros(1,length(PRN));
        Eltc = zeros(1,length(PRN));
        Eer = zeros(1,length(PRN));
        for o = 1:length(PRN)
            d_UDRE(o) = sqrt(I(:,o)'*Cov{1,PRN(o)}*I(:,o)) + Ec(o);
        end
        d_UDRE(isnan(d_UDRE)) = 1;
        Sig2_flt = (sqrt(Sig2_UDRE(PRN)).*d_UDRE + Efc + Errc + Eltc + Eer).^2; % (A.4.5.1 of RTCA DO-229D)
    else
        Sig2_flt = Sig2_UDRE(PRN); % (J.2.2 of RTCA DO-229D)
    end
    Sig2_air = 0.15^2 + (0.13 + 0.53*exp(-Ele/10)).^2; % Airborn equipment (J.2.4 of RTCA DO-229D)
    Sig2_tro = (0.12.*(1.001./sqrt(0.002001+sind(Ele).^2))).^2; % Tropospheric vertical error = 0.12 meters (A.4.2.5 of RTCA DO-229D)
    Sig2 = Sig2_flt + Sig2_UIRE + Sig2_air + Sig2_tro ; % Total Sigma^2
    W = diag(1./Sig2); % Weighting matrix
    % Filter with elevations
    Indx_Cut = Ele<Ele_Cut;
    zXYZs(:,Indx_Cut) = [];
    R(Indx_Cut) = [];
    Sig2(Indx_Cut) = [];
    zW = diag(1./Sig2);
else
    Sig2 = ones(1,length(PRN));
    W = diag(1./Sig2); % Weighting matrix
    zW = W;
end

if (length(R) > 3) % Check number of satellite
    cIter = 0;      % Interative count
    nIter = 100;    % Iterative number
    dX = inf;       % Initial dX
    Bias = 0;       % Initial bias
    while (norm(dX) > 10^-3) && (cIter < nIter) % Iterative method
        cIter = cIter + 1;
        
        % Computed line-of-sight vectors and ranges from satellite to receiver
        pj = zXYZs - XYZu(:)*ones(1,size(zXYZs,2));
        Norm_pj = sqrt(sum(pj.^2,1));
        abc = -pj./(ones(3,1)*Norm_pj); % line of sight unit vectors
        
        % Computed the a-priori range
        P_hat = Norm_pj' + Bias;
        P = R' - P_hat;                  % matrix P
        
        A = [abc',ones(length(abc),1)]; % matrix A
        
        % MMSE solution
        dX = (A'*zW*A)\(A'*zW*P);
        
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

    % === Computed protection levels for EN route though Localizer Performance with Vertical (LPV) approch (J.1 of RTCA DO-229D)
    G = [-cosd(Ele).*sind(Azi); -cosd(Ele).*cosd(Azi); -sind(Ele); ones(1,size(Ele,2))]'; % The geometrix G
    if K > 2
        D = inv(G'*W*G);
        Sig2_E = D(1,1);
        Sig2_N = D(2,2);
        Sig_EN = D(2,1);
        Sig_U = (D(3,3))^0.5;
        SD.Mejor = sqrt(((Sig2_E + Sig2_N)/2) + sqrt((((Sig2_E - Sig2_N)/2)^2) + Sig_EN^2));
        SD.U = Sig_U;
    else
        SD = [];
    end
else
    Pos = [];
    SD = [];
end
