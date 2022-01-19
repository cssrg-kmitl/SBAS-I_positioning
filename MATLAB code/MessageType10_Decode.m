function [Output] = MessageType10_Decode(Type10)
%% ==============================================================
% Objective: To convert the message types 10 to the degradation parameters (A.4.4.6 of RTCA DO-229D).
% Example: [Output] = MessageType10_Decode(Type10).
% Type10 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================================
Char_hex = char(Type10)';   % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Output.Brre = bin2dec(Temp(15:24)).*0.002;        % 10-Bits (m)
Output.Cltc_1sb = bin2dec(Temp(25:34)).*0.002;    % 10-Bits (m)
Output.Cltc_v1 = bin2dec(Temp(35:44)).*0.00005;	% 10-Bits (m/s)
Output.Iltc_v1 = bin2dec(Temp(45:53));     % 9-Bits (s)
Output.Cltc_v0 = bin2dec(Temp(54:63)).*0.002;	% 10-Bits (m)
Output.Iltc_v0 = bin2dec(Temp(64:72));     % 9-Bits (s)
Output.Cgeo_1sb = bin2dec(Temp(73:82)).*0.0005;    % 10-Bits (m)
Output.Cgeo_v = bin2dec(Temp(83:92)).*0.00005;    % 10-Bits (m/s)
Output.Igeo = bin2dec(Temp(93:101));     % 9-Bits (s)
Output.Cer = bin2dec(Temp(102:107)).*0.5;     % 6-Bits (m)
Output.Ciono_step = bin2dec(Temp(108:117)).*0.001;    % 10-Bits (m)
Output.Iiono = bin2dec(Temp(118:126));     % 9-Bits (s)
Output.Ciono_ramp = bin2dec(Temp(127:136)).*0.000005;    % 10-Bits (m/s)
Output.RSS_UDRE = bin2dec(Temp(137)); % 1-Bits
Output.RSS_iono = bin2dec(Temp(138)); % 1-Bits
Output.C_cov = bin2dec(Temp(139:145)).*0.1;     % 7-Bits
% Spare = bin2dec(Bin(:,146:226));     % 81-Bits (No used)
clear Temp
