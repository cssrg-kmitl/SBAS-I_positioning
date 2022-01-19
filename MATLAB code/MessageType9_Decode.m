function [Output] = MessageType9_Decode(Type9)
%% ==============================================================
% Objective: To convert the message types 9 to the GEO navigation message (A.4.4.11 of RTCA DO-229D).
% Example: [Output] = MessageType9_Decode(Type9).
% Type9 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================================
Char_hex = char(Type9)';    % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
URA_v = [2 2.8 4 5.7 8 11.3 16 32 64 128 256 512 1024 2048 4096 NaN]; %URA value
Output.Res = bin2dec(Temp(15:22));	% 8-Bits(Reserved)
Output.t_0 = bin2dec(Temp(23:35)).*16;	% 13-Bits: Time of applicability of the message (seconds)
Output.URA = URA_v(bin2dec(Temp(36:39))+1);	% 4-Bits: User Range Accurracy (m)
Bin2Dec = bin2dec(Temp(40:69));	% 30-Bits: Earth Centered Earth Fixed (m)
Output.Xg = (Bin2Dec > 536870911).*(Bin2Dec*0.08 - 0.08*2^30) + (Bin2Dec <= 536870911).*(Bin2Dec*0.08);
Bin2Dec = bin2dec(Temp(70:99));	% 30-Bits: Earth Centered Earth Fixed (m)
Output.Yg = (Bin2Dec > 536870911).*(Bin2Dec*0.08 - 0.08*2^30) + (Bin2Dec <= 536870911).*(Bin2Dec*0.08);
Bin2Dec = bin2dec(Temp(100:124));   % 25-Bits: Earth Centered Earth Fixed (m)
Output.Zg = (Bin2Dec > 16777215).*(Bin2Dec*0.4 - 0.4*2^25) + (Bin2Dec <= 16777215).*(Bin2Dec*0.4);
Bin2Dec = bin2dec(Temp(125:141));	% 17-Bits: Rate of Change (m/s)
Output.Xroc = (Bin2Dec > 65535).*(Bin2Dec*0.000625 - 0.000625*2^17) + (Bin2Dec <= 65535).*(Bin2Dec*0.000625);
Bin2Dec = bin2dec(Temp(142:158));	% 17-Bits: Rate of Change (m/s)
Output.Yroc = (Bin2Dec > 65535).*(Bin2Dec*0.000625 - 0.000625*2^17) + (Bin2Dec <= 65535).*(Bin2Dec*0.000625);
Bin2Dec = bin2dec(Temp(159:176));	% 18-Bits: Rate of Change (m/s)
Output.Zroc = (Bin2Dec > 131071).*(Bin2Dec*0.004 - 0.004*2^18) + (Bin2Dec <= 131071).*(Bin2Dec*0.004);
Bin2Dec = bin2dec(Temp(177:186));	% 10-Bits: Accelaration (m/s^2)
Output.Xacc = (Bin2Dec > 511).*(Bin2Dec*0.0000125 - 0.0000125*2^10) + (Bin2Dec <= 511).*(Bin2Dec*0.0000125);
Bin2Dec = bin2dec(Temp(187:196));	% 10-Bits: Accelaration (m/s^2)
Output.Yacc = (Bin2Dec > 511).*(Bin2Dec*0.0000125 - 0.0000125*2^10) + (Bin2Dec <= 511).*(Bin2Dec*0.0000125);
Bin2Dec = bin2dec(Temp(197:206));	% 10-Bits: Accelaration (m/s^2)
Output.Zacc = (Bin2Dec > 511).*(Bin2Dec*0.0000625 - 0.0000625*2^10) + (Bin2Dec <= 511).*(Bin2Dec*0.0000625);
Bin2Dec = bin2dec(Temp(207:218));	% 12-Bits: Time offset (Seconds)
Output.a_gf0 = (Bin2Dec > 2047).*(Bin2Dec*(2^-31) - (2^-31)*2^12) + (Bin2Dec <= 2047).*(Bin2Dec*(2^-31));
Bin2Dec = bin2dec(Temp(219:226));	% 8-Bits: Time drift (Seconds/s)
Output.a_gf1 = (Bin2Dec > 127).*(Bin2Dec*(2^-40) - (2^-40)*2^8) + (Bin2Dec <= 127).*(Bin2Dec*(2^-40));
clear Temp;
