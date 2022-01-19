%-------------------------------------------------------------------------------------------------------------
% Objective: To simulate the performance of the SBAS-I systems.
% Input: The observation file, navigation file, and SBAS file of the RINEX version 2 format.
% Output: Receiver position and standard deviations of the "Mejor" and "U" to compute the HPL and VPL.
% Experimental setup:   1) Fill the RINEX file name.
%                       2) Chosen one of the SBAS systems within the "List" variable.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (April 2020).
%-------------------------------------------------------------------------------------------------------------
clear; % Claer memory
tic; % Count time
%% === Reading RINEX files ===
File_O = 'KMIT0600.19o'; % Observation file name (String type)
File_N = 'KMIT0600.19n'; % Navigation file name (String type)
File_S = 'KMIT0600.19s'; % SBAS file name (String type)
[Obs.Date.St,Obs.Ep.G,Obs.Type.G,Obs.Unit.G,Obs.PRN.G,Obs.rcv,Obs.Data.G,Obs.PRN.G,Obs.rcvpos,Obs.antdel,Obs.anttype,Obs.Com] =...
    readrinexobs(File_O);
[Nav.Sats,Nav.rcvname,Nav.Eph.G,Nav.PRN.G,Nav.ion,Nav.dutc,Nav.com] = readrinexnav(File_N);
Nav.Eph.G = Nav.Eph.G';
Nav.Ion.GPSA = Nav.ion(1,:);
Nav.Ion.GPSB = Nav.ion(2,:);
% === Check date
Date = datetime(Obs.Date.St(1),Obs.Date.St(2),Obs.Date.St(3));
DOW = day(Date,'dayofweek');    % Day of week
DOY = day(Date,'dayofyear');    % Day of year
t0 = (DOW-1)*86400;             % Start-time on a day
% === Read SBAS information
[List] = SpliteSBASsys(File_S);
SBAS_Data = List.GAGAN_Data; % **** Chosen one of the SBAS systems ****
SBAS_PRN = List.GAGAN_PRN;   % **** Chosen one of the SBAS systems ****
for n = 1:length(SBAS_PRN)
    L_Ep(n) = length(SBAS_Data{n, 1});
end
[~,Indx_PRN] = max(L_Ep); % Chosen the most data size

%% === Parameter-setup ===
GPSconstant;                % Loaded constant parameters
Sampling = 6;               % Sampling rate (Integer type)
End_time = 86399;           % End time of a day
Time = 0:Sampling:End_time; % Computed times
Epoch = 0:End_time;         % Created the new Epoch to sample the SBAS corrections
L_Time = length(Time);      % Length of time matrix

%% === SPS output variables ===
SPS.LLA = nan(L_Time,3);
SPS.XYZ = nan(L_Time,3);
SPS.Sat_Num = nan(L_Time,1);
SODs = nan(L_Time,1);

%% === SPS+SBAS output variables ===
Row = 1;
Col = 51; % PRN number aided by SBAS system
% === Fast&Long-term corrections
FL.LLA = nan(L_Time,3);
FL.XYZ = nan(L_Time,3);
FL.Sat_Num = nan(L_Time,1);
FL.SD_Mejor = nan(L_Time,1);
FL.SD_U = nan(L_Time,1);
% === Fast&Long-term corrections with Ionospheric grid point (IGP)
FLI.LLA = nan(L_Time,3);
FLI.XYZ = nan(L_Time,3);
FLI.Sat_Num = nan(L_Time,1);
FLI.SD_Mejor = nan(L_Time,1);
FLI.SD_U = nan(L_Time,1);

%% === Buffer variables of SBAS parameters ===
IODP = nan(Row,8);              % Issue of data PRN mask
IODF = nan(Row,4);              % Issue of data fast corrections
IODI = nan(Row,2);              % Issue of data IGP
IOD_Buf = nan(Row,Col);         % IOD buffer
% === MT1 (A.4.4.2 of RTCA DO-229D)
PRN_Mask = nan(Row,Col);        % PRN mask number
% === MT2-5 (A.4.4.3 of RTCA DO-229D)
PRC = nan(Row,Col);             % Fast corrections
oPRC = nan(Row,Col);            % Previous fast corrections
UDRE = nan(Row,Col);            % User differential range errors
Sig2_UDRE = nan(Row,Col);       % Variance of UDREs
% === MT7 (A.4.4.5 of RTCA DO-229D)
ai2 = ones(Row,Col)*0;          % Fast correction degradation factors (m/s^2)
Ifc1 = ones(Row,Col)*180;       % Time-out interval (En-route through LNAV approach)
Ifc2 = ones(Row,Col)*120;       % Time-out interval (LNAV/VNAV,LPV,LP approach)
t_Update = ones(Row,Col)*60;    % Fast corrections updated interval (sec)
t_lat = nan(Row,Col);           % System latency (sec)
% === MT9 (A.4.4.11 of RTCA DO-229D)
t_0 = nan(Row,1);               % Time of applicability
URA = nan(Row,1);               % Accuracy exponent
XYZ_G = nan(Row,3);             % ECEF with WGS-84 of the geostationary satellite
rXYZ_G = nan(Row,3);            % Rate-of-change of the XYZ_G
aXYZ_G = nan(Row,3);            % Accelleration of the XYZ_G
a_Gf0 = nan(Row,1);             % Time offset
a_Gf1 = nan(Row,1);             % Time drift
% === MT10 (A.4.4.6 of RTCA DO-229D)
% Degradation factors
Brre = nan(Row,1);
Cltc_1sb = nan(Row,1);
Cltc_v1 = nan(Row,1);
Iltc_v1 = nan(Row,1);
Cltc_v0 = nan(Row,1);
Iltc_v0 = nan(Row,1);
Cgeo_1sb = nan(Row,1);
Cgeo_v = nan(Row,1);
Igeo = nan(Row,1);
Cer = nan(Row,1);
Ciono_step = nan(Row,1);
Iiono = nan(Row,1);
Ciono_ramp = nan(Row,1);
RSSu = nan(Row,1);
RSSi = nan(Row,1);
C_cov = nan(Row,1);
% === MT17 (A.4.4.12 of RTCA DO-229D)
% Geo almanacs message parameters
DataID = nan(Row,3);
PRN_Geo = nan(Row,3);       % PRN number of the geostationary satellites
ID = nan(Row,3);            % Service provider ID
Res = nan(Row,3);           % Reserved
Bro = nan(Row,3);           % Broadcast
Cor = nan(Row,3);           % Correction
Ran = nan(Row,3);           % Ranging
XYZ_Geo = cell(Row,1);      % Rate-of-change of the XYZ_G
rXYZ_Geo = cell(Row,1);     % Rate-of-change of the XYZ_G
% === MT25 (A.4.4.7 of RTCA DO-229D)
Vcod = nan(Row,Col);        % Velocity codes
d_XYZ = cell(Row,Col);      % Errors of the ECEF with WGS-84
d_XYZr = cell(Row,Col);     % Rate-of-change of d_XYZ
D_af0 = nan(Row,Col);       % Time offset error corrections
D_af1 = nan(Row,Col);       % Time drift error corrections
t_o = nan(Row,Col);         % Time-of-day applicability
% === MT18 (A.4.4.9 of RTCA DO-229D)
BandLat = cell(Row,11);     % Band of the latitude axist 
BandLon = cell(Row,11);     % Band of the longitude axist
% === MT26 (A.4.4.10 of RTCA DO-229D)
BandIVDE = cell(Row,11);    % Bands of the IGP Vertical Delay Estimate
IVDE = nan(Row,201,11);     % Buffer of the IVDE values
BandGIVE = cell(Row,11);    % Bands of the Grid Ionospheric Vertical Errors
GIVE = nan(Row,201,11);     % Buffer of the GIVE values
BandSig2_GIVE = cell(Row,11);   % Bands of the GIVE variance
Sig2_GIVE = nan(Row,201,11);	% Buffer of the GIVE variance values
% === MT28 (A.4.4.16 of RTCA DO-229D)
Cov = cell(Row,Col);        % Covariance matrix
SF = nan(Row,Col);          % Scale factor
% Initail values for the cell data types
k = Row;
XYZ_Geo{k} = nan(3,3);
rXYZ_Geo{k} = nan(3,3);
for o = 1:Col
    d_XYZ{k,o} = nan(3,1);
    d_XYZr{k,o} = nan(3,1);
    Cov{k,o} = nan(4);
end
for o = 1:11
    BandLat{k,o} = nan(1,201);
    BandLon{k,o} = nan(1,201);
    BandIVDE{k,o} = nan(1,201);
    BandGIVE{k,o} = nan(1,201);
    BandSig2_GIVE{k,o} = nan(1,201);
end

%% === Logged variables of SBAS parameters ===
SBAS.PRN = [];
SBAS.PRC = [];
SBAS.oPRC = [];
SBAS.UDRE = [];
SBAS.Sig2_UDRE = [];
SBAS.Time_pr = [];
SBAS.Time_ap = [];
SBAS.t_lat = [];
SBAS.ai2 = [];
SBAS.Ifc1 = [];
SBAS.Ifc2 = [];
SBAS.t_Update = [];
SBAS.Brre = [];
SBAS.Cltc_1sb = [];
SBAS.Cltc_v1 = [];
SBAS.Iltc_v1 = [];
SBAS.Cltc_v0 = [];
SBAS.Iltc_v0 = [];
SBAS.Cgeo_1sb = [];
SBAS.Cgeo_v = [];
SBAS.Igeo = [];
SBAS.Cer = [];
SBAS.Ciono_step = [];
SBAS.Iiono = [];
SBAS.Ciono_ramp = [];
SBAS.RSSu = [];
SBAS.RSSi = [];
SBAS.C_cov = [];
SBAS.Vcod = [];
SBAS.d_XYZ = [];
SBAS.d_XYZr = [];
SBAS.d_af0 = [];
SBAS.d_af1 = [];
SBAS.t_o = [];
SBAS.IOD = [];
SBAS.Lat = [];
SBAS.Lon = [];
SBAS.IVDE = [];
SBAS.GIVE = [];
SBAS.Sig2_GIVE = [];
SBAS.Cov = [];
SBAS.SF = [];

%% === Smoothed variables ===
alph = Sampling/100;
LambdaL1 = GPScons.lambda1;
L1_Pre = zeros(1,Col); % Previous pseoranges of carrier phase
Pn_Pre = zeros(1,Col); % Previous pseoranges of code

if (str2double(SBAS_Data{Indx_PRN,1}(1,2))- t0) == Time(1) % Check the matched data with the start time of a day
    %% === Initail values ===
    Height_SPS = 0;         % Estimated user altitude based on the Standard Positioning System (SPS)
    XYZu_SPS = [0; 0; 0];   % Estimated user position based on the SPS
    Height_FL = 0;          % Estimated user altitude based on the SPS + SBAS_FL
    XYZu_FL = [0; 0; 0];    % Estimated User position based on the SPS + SBAS_FL
    Height_FLI = 0;         % Estimated user altitude based on the SPS + SBAS_FLI
    XYZu_FLI = [0; 0; 0];   % Estimated User position based on the SPS + SBAS_FLI
    G1 = 1;     % Computed counter of SPS
    S1 = 1;     % Computed counter of SPS+SBAS_FL
    S2 = 1;     % Computed counter of SPS+SBAS_FLI
    C = 0;      % Availability check flag of the SBAS parameters
    T = 1;      % Index counter
    k = 1;      % Row
    IODP(1) = 4; % IODP flag
    IODI(1) = 4; % IODI flag
    Time_out = zeros(Row,Col);  % Time-out buffer
    Time_pr = zeros(1,Col);     % Previous time buffer
    Time_ap = ones(1,Col);      % Applicability time buffer
    
    for t = 1:End_time
        if (t<L_Ep(Indx_PRN))
            Indx_Ep = find(Epoch == str2double(SBAS_Data{Indx_PRN,1}(t,2))-t0);
            if (Indx_Ep<=L_Ep(Indx_PRN))
                Ch = str2double(SBAS_Data{Indx_PRN,1}(Indx_Ep,4));
                %% === Observed the SBAS parameters ===
                if (Ch==1)
                    [MT1] = MessageType1_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(1) = MT1.IODP;
                    PRN_Mask(k,MT1.PRN_Mask) = MT1.PRN_Mask;
                elseif (Ch==2)
                    [MT2] = MessageType2to5_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(2) = MT2.IODP;
                    IODF(1) = MT2.IODF;
                    oPRC(k,1:13) = PRC(k,1:13);
                    PRC(k,1:13) = MT2.PRC;
                    UDRE(k,1:13) = MT2.UDRE;
                    Sig2_UDRE(k,1:13) = MT2.Sig2_UDRE;
                    Time_out(k,1:13)= 0;
                    Time_pr(1:13) = Time_ap(1:13);
                    Time_ap(1:13) = t;
                elseif (Ch==3)
                    [MT3] = MessageType2to5_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(3) = MT3.IODP;
                    IODF(2) = MT3.IODF;
                    oPRC(k,14:26) = PRC(k,14:26);
                    PRC(k,14:26) = MT3.PRC;
                    UDRE(k,14:26) = MT3.UDRE;
                    Sig2_UDRE(k,14:26) = MT3.Sig2_UDRE;
                    Time_out(k,14:26) = 0;
                    Time_pr(14:26) = Time_ap(14:26);
                    Time_ap(14:26) = t;
                elseif (Ch==4)
                    [MT4] = MessageType2to5_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(4) = MT4.IODP;
                    IODF(3) = MT4.IODF;
                    oPRC(k,27:39) = PRC(k,27:39);
                    PRC(k,27:39) = MT4.PRC;
                    UDRE(k,27:39) = MT4.UDRE;
                    Sig2_UDRE(k,27:39) = MT4.Sig2_UDRE;
                    Time_out(k,27:39) = 0;
                    Time_pr(27:39) = Time_ap(27:39);
                    Time_ap(27:39) = t;
                elseif (Ch==7)
                    [MT7] = MessageType7_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(5) = MT7.IODP;
                    t_lat(k,:) = MT7.t_lat;
                    ai2(k,:) = MT7.ai2;
                    Ifc1(k,:) = MT7.Ifc1;
                    Ifc2(k,:) = MT7.Ifc2;
                    t_Update(k,:) = MT7.t_Update;
                elseif (Ch==9)
                    [MT9] = MessageType9_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    t_0(k) = MT9.t_0;
                    URA(k) = MT9.URA;
                    XYZ_G(k,:) = [MT9.Xg MT9.Yg MT9.Zg];
                    rXYZ_G(k,:) = [MT9.Xroc MT9.Yroc MT9.Zroc];
                    aXYZ_G(k,:) = [MT9.Xacc MT9.Yacc MT9.Zacc];
                    a_Gf0(k) = MT9.a_gf0;
                    a_Gf1(k) = MT9.a_gf1;
                elseif (Ch==10)
                    [MT10] = MessageType10_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    Brre(k) = MT10.Brre;
                    Cltc_1sb(k) = MT10.Cltc_1sb;
                    Cltc_v1(k) = MT10.Cltc_v1;
                    Iltc_v1(k) = MT10.Iltc_v1;
                    Cltc_v0(k) = MT10.Cltc_v0;
                    Iltc_v0(k) = MT10.Iltc_v0;
                    Cgeo_1sb(k) = MT10.Cgeo_1sb;
                    Cgeo_v(k) = MT10.Cgeo_v;
                    Igeo(k) = MT10.Igeo;
                    Cer(k) = MT10.Cer;
                    Ciono_step(k) = MT10.Ciono_step;
                    Iiono(k) = MT10.Iiono;
                    Ciono_ramp(k) = MT10.Ciono_ramp;
                    RSSu(k) = MT10.RSS_UDRE;
                    RSSi(k) = MT10.RSS_iono;
                    C_cov(k) = MT10.C_cov;
                elseif (Ch==17)
                    [MT17] = MessageType17_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    DataID(k,:) = MT17.DataID;
                    PRN_Geo(k,:) = MT17.PRN_num;
                    ID(k,:) = MT17.Health.ID;
                    Res(k,:) = MT17.Health.Res;
                    Bro(k,:) = MT17.Health.Bro;
                    Cor(k,:) = MT17.Health.Cor;
                    Ran(k,:) = MT17.Health.Ran;
                    XYZ_Geo{k,1} = [MT17.Xg; MT17.Yg; MT17.Zg]';
                    rXYZ_Geo{k,1} = [MT17.Xroc; MT17.Yroc; MT17.Zroc]';
                elseif (Ch==18)
                    [MT18] = MessageType18_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODI(1) = MT18.IODI;
                    BandLat{k,MT18.Band+1} = MT18.Lat;
                    BandLon{k,MT18.Band+1} = MT18.Lon;
                elseif (Ch==25)
                    [MT25] = MessageType25_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODP(6) = MT25.IODP(1);
                    IODP(7) = MT25.IODP(2);
                    Vc = MT25.Vcod;
                    PRN_mask = MT25.PRN_Mask;
                    IOD = MT25.IOD;
                    dX = MT25.dX; dY = MT25.dY; dZ = MT25.dZ;
                    dXr = MT25.dXr; dYr = MT25.dYr; dZr = MT25.dZr;
                    d_af0 = MT25.d_af0; d_af1 = MT25.d_af1;
                    to = MT25.t_o;
                    PRN_mask(PRN_mask==0) = [];
                    for o = 1:length(PRN_mask)
                        d_XYZ{k,PRN_mask(o)} = [dX(o);dY(o);dZ(o)];
                        d_XYZr{k,PRN_mask(o)} = [dXr(o);dYr(o);dZr(o)];
                    end
                    Vcod(k,PRN_mask) = Vc(PRN_mask>0);
                    IOD_Buf(k,PRN_mask) = IOD(PRN_mask>0);
                    D_af0(k,PRN_mask) = d_af0(PRN_mask>0);
                    D_af1(k,PRN_mask) = d_af1(PRN_mask>0);
                    t_o(k,PRN_mask) = to(PRN_mask>0);
                elseif (Ch==26)
                    [MT26] = MessageType26_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    IODI(2) = MT26.IODI;
                    Band_numb = MT26.Band_numb;
                    Block_ID = MT26.Block_ID;
                    buffIGP_VDE = MT26.IGP_VDE;
                    buffGIVE = MT26.GIVE;
                    bSig2_GIVE = MT26.Sig2_GIVE;
                    if sum(~isnan(BandLat{k,Band_numb+1}))>0
                        nIGP_Block = find(~isnan(BandLat{k,Band_numb+1}));
                        L26 = length(buffIGP_VDE);
                        Temp = 1+Block_ID*15:L26+Block_ID*15;
                        Temp(Temp>length(nIGP_Block)) = []; % Cutting the tail array is over
                        Index = nIGP_Block(Temp);
                        IVDE(k,Index,Band_numb+1) = buffIGP_VDE(1:length(Index))';
                        GIVE(k,Index,Band_numb+1) = buffGIVE(1:length(Index))';
                        Sig2_GIVE(k,Index,Band_numb+1) = bSig2_GIVE(1:length(Index))';
                        BandIVDE{k,Band_numb+1} = IVDE(k,:,Band_numb+1);
                        BandGIVE{k,Band_numb+1} = GIVE(k,:,Band_numb+1);
                        BandSig2_GIVE{k,Band_numb+1} = Sig2_GIVE(k,:,Band_numb+1);
                    end
                elseif (Ch==28)
                    [MT28] = MessageType28_Decode(SBAS_Data{Indx_PRN,1}(Indx_Ep,5));
                    if IODP(1) == MT28.IODP
                        IODP(8) = MT28.IODP;
                        if (MT28.PRN1>0)
                            Cov{k,MT28.PRN1} = MT28.Cov1;
                            SF(k,MT28.PRN1) = MT28.SF1;
                        end
                        if (MT28.PRN2>0)
                            Cov{k,MT28.PRN2} = MT28.Cov2;
                            SF(k,MT28.PRN2) = MT28.SF2;
                        end
                    end
                end
            end
        end
        
        %% === Logged the SBAS parameters ===
        if sum(IODP(:)==IODP(1))>7
            SBAS.PRN = PRN_Mask(k,:);   % PRN mask numbers
            SBAS.Time_pr = Time_pr;     % Time of applicability of a previous fast corrections
            SBAS.Time_ap = Time_ap;     % Time of applicability of the most recent fast corrections
            SBAS.oPRC = oPRC(k,:);      % Previous PRC
            if sum(IODF(:)<3)>2
                SBAS.PRC = PRC(k,:);
                SBAS.UDRE = UDRE(k,:);
                SBAS.Sig2_UDRE = Sig2_UDRE(k,:);
            end
            SBAS.t_lat = t_lat(k,:);
            SBAS.ai2 = ai2(k,:);
            SBAS.Ifc1 = Ifc1(k,:);
            SBAS.Ifc2 = Ifc2(k,:);
            SBAS.t_Update = t_Update(k,:);
            SBAS.Brre = Brre(k);
            SBAS.Cltc_1sb = Cltc_1sb(k);
            SBAS.Cltc_v1 = Cltc_v1(k);
            SBAS.Iltc_v1 = Iltc_v1(k);
            SBAS.Cltc_v0 = Cltc_v0(k);
            SBAS.Iltc_v0 = Iltc_v0(k);
            SBAS.Cgeo_1sb = Cgeo_1sb(k);
            SBAS.Cgeo_v = Cgeo_v(k);
            SBAS.Igeo = Igeo(k);
            SBAS.Cer = Cer(k);
            SBAS.Ciono_step = Ciono_step(k);
            SBAS.Iiono = Iiono(k);
            SBAS.Ciono_ramp = Ciono_ramp(k);
            SBAS.RSSu = RSSu(k);
            SBAS.RSSi = RSSi(k);
            SBAS.C_cov = C_cov(k);
            SBAS.Vcod = Vcod(k,:);
            SBAS.d_XYZ = d_XYZ(k,:);
            SBAS.d_XYZr = d_XYZr(k,:);
            SBAS.d_af0 = D_af0(k,:);
            SBAS.d_af1 = D_af1(k,:);
            SBAS.t_o = t_o(k,:);
            SBAS.IOD = IOD_Buf(k,:);
            SBAS.Cov = Cov(k,:);
            SBAS.SF = SF(k,:);
            C = 1;
            if (sum(IODI(:)== IODI(1))>1)
                SBAS.Lat = BandLat(k,:);
                SBAS.Lon = BandLon(k,:);
                SBAS.IVDE = BandIVDE(k,:);
                SBAS.GIVE = BandGIVE(k,:);
                SBAS.Sig2_GIVE  = BandSig2_GIVE(k,:);
            end
            if isempty(SBAS.Lat)
                C = 0;
            end
        end
        if ~mod(t,Sampling)
            %% === Estimated the receiver position ===
            [Smoothed] = CarrierPhaseSmooth(Obs,4,'C1','L1',Pn_Pre,L1_Pre,LambdaL1,alph,t);
            Pn_Pre = Smoothed.Pn_Pre;
            L1_Pre = Smoothed.L1_Pre;
            PRN = Smoothed.PRN;     % PRN number of the GPS satellites
            SOD = Smoothed.SOD;     % Considered time of each GPS satellite
            Pn = Smoothed.Pn;       % Pseudoranges between antennas of satellite and receiver smoothed by the carrier phase
            if (length(PRN)<4)
                continue;
            end
            
            if C
                % === Computed the receiver position based on the SPS+SBAS_FL
                [Pos_FL] = ComputedPositionByGPSwithSBAS_FL(DOY,SOD,PRN,Obs,Nav,Pn,XYZu_FL,Height_FL,SBAS,GPScons,S1);
                if ~isempty(Pos_FL)
                    XYZu_FL = Pos_FL.xyz';
                    Height_FL = Pos_FL.llh(3);
                    FL.LLA(T,:) = Pos_FL.llh;
                    FL.XYZ(T,:) = Pos_FL.xyz';
                    FL.Sat_Num(T) = Pos_FL.Sat_Num;
                    S1 = S1+1;
                end
                
                % === Computed the receiver position based on the SPS+SBAS_FLI
                [Pos_FLI,SD_FLI] = ComputedPositionByGPSwithSBAS_FLI(DOY,SOD,PRN,Obs,Nav,Pn,XYZu_FLI,Height_FLI,SBAS,GPScons,S2);
                if ~isempty(Pos_FLI)
                    XYZu_FLI = Pos_FLI.xyz';
                    Height_FLI = Pos_FLI.llh(3);
                    FLI.LLA(T,:) = Pos_FLI.llh;
                    FLI.XYZ(T,:) = Pos_FLI.xyz';
                    FLI.Sat_Num(T) = Pos_FLI.Sat_Num;
                    if ~isempty(SD_FLI)
                        FLI.SD_Mejor(T) = SD_FLI.Mejor;
                        FLI.SD_U(T) = SD_FLI.U;
                    end
                    S2 = S2+1;
                end
                
                % === Computed the receiver position based on the SPS
                [Pos_SPS] = ComputedPositionByGPS(DOY,SOD,PRN,Obs,Nav,Pn,XYZu_SPS(:),Height_SPS,GPScons,G1);
                if ~isempty(Pos_SPS)
                    XYZu_SPS = Pos_SPS.xyz';
                    Height_SPS = Pos_SPS.llh(3);
                    SPS.LLA(T,:) = Pos_SPS.llh;
                    SPS.XYZ(T,:) = Pos_SPS.xyz';
                    SPS.Sat_Num(T) = Pos_SPS.Sat_Num;
                    G1 = G1+1;
                end
                SODs(T) = SOD(1);
            end
            T = T+1;
        end
    end
else
    error('Data is mismatch');
end
toc; % Display count time
PlotSimulationResults;
