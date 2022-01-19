function [Output] = MessageType17_Decode(Type17)
%% ==============================================================
% Objective: To convert the message types 17 to the GEO satellite almanacs (A.4.4.12 of RTCA DO-229D).
% Example: [Output] = MessageType17_Decode(Type17).
% Type17 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================================
Char_hex = char(Type17)';   % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Output.DataID = [bin2dec(Temp(15:16)) bin2dec(Temp(82:83)) bin2dec(Temp(149:150))];	% 3 x 2-Bits(Reserved)
Output.PRN_num = [bin2dec(Temp(17:24)) bin2dec(Temp(84:91)) bin2dec(Temp(151:158))];	% 3 x 8-Bits
Health.ID = [bin2dec(Temp(25:28)) bin2dec(Temp(92:95)) bin2dec(Temp(159:162))]; % 3 x 4-Bits: Health.Service provider ID
Health.Res = [bin2dec(Temp(29)) bin2dec(Temp(96)) bin2dec(Temp(163))];   % 3 x 1-Bits: Health.Reserved
Health.Bro = [bin2dec(Temp(30)) bin2dec(Temp(97)) bin2dec(Temp(164))];   % 3 x 1-Bits: Health.Broadcast integrity
Health.Cor = [bin2dec(Temp(31)) bin2dec(Temp(98)) bin2dec(Temp(165))];   % 3 x 1-Bits: Health.Correction
Health.Ran = [bin2dec(Temp(32)) bin2dec(Temp(99)) bin2dec(Temp(166))];	 % 3 x 1-Bits: Health.Ranging
Output.Health = Health;
Bin2Dec = [bin2dec(Temp(33:47)) bin2dec(Temp(100:114)) bin2dec(Temp(167:181))]; % 3 x 15-Bits: Earth Centered Eart Fixed (m)
Output.Xg = (Bin2Dec > 16383).*(Bin2Dec*2600 - 85196800) + (Bin2Dec <= 16383).*(Bin2Dec*2600);
Bin2Dec = [bin2dec(Temp(48:62)) bin2dec(Temp(115:129)) bin2dec(Temp(182:196))];	% 3 x 15-Bits: Earth Centered Eart Fixed (m)
Output.Yg = (Bin2Dec > 16383).*(Bin2Dec*2600 - 85196800) + (Bin2Dec <= 16383).*(Bin2Dec*2600);
Bin2Dec = [bin2dec(Temp(63:71)) bin2dec(Temp(130:138)) bin2dec(Temp(197:205))];	% 3 x 9-Bits:  Earth Centered Eart Fixed (m)
Output.Zg = (Bin2Dec > 255).*(Bin2Dec*26000 - 13312000) + (Bin2Dec <= 255).*(Bin2Dec*26000);
Bin2Dec = [bin2dec(Temp(72:74)) bin2dec(Temp(139:141)) bin2dec(Temp(206:208))];	% 3 x 3-Bits: Rat of Change (m/s)
Output.Xroc = (Bin2Dec > 3).*(Bin2Dec*10 - 80) + (Bin2Dec <= 3).*(Bin2Dec*10);
Bin2Dec = [bin2dec(Temp(75:77)) bin2dec(Temp(142:144)) bin2dec(Temp(209:211))];	% 3 x 3-Bits: Rat of Change (m/s)
Output.Yroc = (Bin2Dec > 3).*(Bin2Dec*10 - 80) + (Bin2Dec <= 3).*(Bin2Dec*10);
Bin2Dec = [bin2dec(Temp(78:81)) bin2dec(Temp(145:148)) bin2dec(Temp(212:215))];	% 3 x 4-Bits: Rat of Change (m/s)
Output.Zroc = (Bin2Dec > 7).*(Bin2Dec*60 - 960) + (Bin2Dec <= 7).*(Bin2Dec*60);
Output.t_o = bin2dec(Temp(216:226))*64; % Tims of day (s)
clear Temp
