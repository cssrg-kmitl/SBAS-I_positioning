function [Output] = MessageType28_Decode(Type28)
%% ====================================================================
% Objective: To convert the message types 28 to the clock-ephemeris covariance matrix (A.4.4.16 of RTCA DO-229D).
% Example: [Output] = MessageType28_Decode(Type28).
% Type28 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ====================================================================
E1 = zeros(4);
E2 = zeros(4);
Char_hex = char(Type28)';	% A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Output.IODP = bin2dec(Temp(15:16));	% IODP status for 2-bits
% ==========================================
Output.PRN1 = bin2dec(Temp(17:22));	% PRN mask for 6-bits
Scale1 = bin2dec(Temp(23:25));	% Scale exponent for 3-bits
% === Diagonal matrix elements of the E1 ===
E1(1,1) = bin2dec(Temp(26:34));	% 9-bits
E1(2,2) = bin2dec(Temp(35:43));	% 9-bits
E1(3,3) = bin2dec(Temp(44:52));	% 9-bits
E1(4,4) = bin2dec(Temp(53:61));	% 9-bits
% === Upper & lower diagonal matrix elements of the E1 ====
Bin2Dec = bin2dec(Temp(62:71));	% 10-bits
E1(1,2) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(72:81));     % 10-bits
E1(1,3) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(82:91));     % 10-bits
E1(1,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(92:101));    % 10-bits
E1(2,3) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(102:111));   % 10-bits
E1(2,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(112:121));   % 10-bits
E1(3,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
% =========================================================
Output.PRN2 = bin2dec(Temp(122:127));	% PRN mask for 6-bits
Scale2 = bin2dec(Temp(128:130));	% Scale exponent for 3-bits
% === Diagonal matrix elements of the E2 ===
E2(1,1) = bin2dec(Temp(131:139));	% 9-bits
E2(2,2) = bin2dec(Temp(140:148));	% 9-bits
E2(3,3) = bin2dec(Temp(149:157));   % 9-bits
E2(4,4) = bin2dec(Temp(158:166));   % 9-bits
% === Upper & lower diagonal matrix elements of the E2 ====
Bin2Dec = bin2dec(Temp(167:176));   % 10-bits
E2(1,2) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(177:186));   % 10-bits
E2(1,3) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(187:196));   % 10-bits
E2(1,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(197:206));   % 10-bits
E2(2,3) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec =  bin2dec(Temp(207:216));  % 10-bits
E2(2,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
Bin2Dec = bin2dec(Temp(217:226));   % 10-bits
E2(3,4) = (Bin2Dec > 511).*(Bin2Dec - 1024) + (Bin2Dec <= 511).*(Bin2Dec);
% =========================================================
clear Temp
SF1 = 2^(Scale1-5); % Scale factor for the E1
SF2 = 2^(Scale2-5); % Scale factor for the E2
R1 = E1*SF1;
R2 = E2*SF2;
Output.SF1 = SF1;
Output.SF2 = SF2;
Output.Cov1 = R1'*R1;	% Covariance matrix for the E1
Output.Cov2 = R2'*R2;	% Covariance matrix for the E2
