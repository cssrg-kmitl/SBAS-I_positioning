function [Output] = MessageType25_Decode(Type25)
%% =======================================================================
% Objective: To convert the message types 25 to the long-term corrections (A.4.4.7 of RTCA DO-229D).
% Example: [Output] = MessageType25_Decode(Type25).
% Type25 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% =======================================================================
% === Parameter variables ===
PRN_Mask = []; % PRN mask no.
IOD = [];   % Issue of data
dX = [];    % ECEF
dY = [];    % ECEF
dZ = [];    % ECEF
dXr = [];   % Rat of change for ECEF
dYr = [];   % Rat of change for ECEF
dZr = [];   % Rat of change for ECEF
d_af0 = []; % Clock offset error corrections
d_af1 = []; % Clock drift error corrections
t_o = [];    % Tim of day availability
IODP = [];  % IODP status
% Spare = [];   
% ===========================
Char_hex = char(Type25)';   % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Vcod = [bin2dec(Temp(:,15)) bin2dec(Temp(:,121))];	% Velocity code of 1-Bit
Output.Vcod = Vcod;
for i = 1:2
    if Vcod(i) == 1 % Velocity code condition is 1
        PRN_Mask = [PRN_Mask bin2dec(Temp((16:21)+(i-1)*106))];    % PRN mask number for 6-Bits
        IOD = [IOD bin2dec(Temp((22:29)+(i-1)*106))];	% 8-Bits
        Bin2Dec = bin2dec(Temp((30:40)+(i-1)*106));
        Real = (Bin2Dec > 1023).*(Bin2Dec*0.125 - 256) + (Bin2Dec <= 1023).*(Bin2Dec*0.125);
        dX = [dX Real]; % 11(m,m/s)
        Bin2Dec = bin2dec(Temp((41:51)+(i-1)*106));
        Real = (Bin2Dec > 1023).*(Bin2Dec*0.125 - 256) + (Bin2Dec <= 1023).*(Bin2Dec*0.125);
        dY = [dY Real];	% 11(m,m/s)
        Bin2Dec = bin2dec(Temp((52:62)+(i-1)*106));
        Real = (Bin2Dec > 1023).*(Bin2Dec*0.125 - 256) + (Bin2Dec <= 1023).*(Bin2Dec*0.125);
        dZ = [dZ Real]; % 11(m,m/s)
        Bin2Dec = bin2dec(Temp((74:81)+(i-1)*106));
        Real = (Bin2Dec > 127).*(Bin2Dec*(2^-11) - 0.125) + (Bin2Dec <= 127).*(Bin2Dec*2^-11);
        dXr = [dXr Real]; % 8-Bits (m,m/s)
        Bin2Dec = bin2dec(Temp((82:89)+(i-1)*106));
        Real = (Bin2Dec > 127).*(Bin2Dec*(2^-11) - 0.125) + (Bin2Dec <= 127).*(Bin2Dec*2^-11);
        dYr = [dYr Real]; % 8-Bits (m,m/s)
        Bin2Dec = bin2dec(Temp((90:97)+(i-1)*106));
        Real = (Bin2Dec > 127).*(Bin2Dec.*(2^-11) - 0.125) + (Bin2Dec <= 127).*(Bin2Dec.*2^-11);
        dZr = [dZr Real]; % 8-Bits (m,m/s)
        Bin2Dec = bin2dec(Temp((63:73)+(i-1)*106));
        Real = (Bin2Dec > 1023).*(Bin2Dec*(2^-31) - (2^-31)*2^11) + (Bin2Dec <= 1023).*(Bin2Dec*2^-31);
        d_af0 = [d_af0 Real]; % 11 (second,second/s)
        Bin2Dec = bin2dec(Temp((98:105)+(i-1)*106));
        Real = (Bin2Dec > 127).*(Bin2Dec*(2^-39) - (2^-39)*2^8) + (Bin2Dec <= 127).*(Bin2Dec*2^-39);
        d_af1 = [d_af1 Real]; % 8-Bits (second,second/s)
        t_o = [t_o bin2dec(Temp((106:118)+(i-1)*106))*16];     % 13-Bits (s)
        IODP = [IODP bin2dec(Temp((119:120)+(i-1)*106))];    % 2-Bits
    else	% Velocity code condition is 0
        PRN_Mask = [PRN_Mask bin2dec(Temp((16:21)+(i-1)*106)) bin2dec(Temp((67:72)+(i-1)*106))];    % 6,6-Bits
        IOD = [IOD bin2dec(Temp((22:29)+(i-1)*106)) bin2dec(Temp((73:80)+(i-1)*106))];	% 8,8-Bits
        Bin2Dec = [bin2dec(Temp((30:38)+(i-1)*106)) bin2dec(Temp((81:89)+(i-1)*106))];
        Real = (Bin2Dec > 255).*(Bin2Dec.*0.125 - 64) + (Bin2Dec <= 255).*(Bin2Dec.*0.125);
        dX = [dX Real];   % 9,9-Bits (m)
        Bin2Dec = [bin2dec(Temp((39:47)+(i-1)*106)) bin2dec(Temp((90:98)+(i-1)*106))];
        Real = (Bin2Dec > 255).*(Bin2Dec.*0.125 - 64) + (Bin2Dec <= 255).*(Bin2Dec.*0.125);
        dY = [dY Real];	% 9,9-Bits (m)
        Bin2Dec = [bin2dec(Temp((48:56)+(i-1)*106)) bin2dec(Temp((99:107)+(i-1)*106))];
        Real = (Bin2Dec > 255).*(Bin2Dec.*0.125 - 64) + (Bin2Dec <= 255).*(Bin2Dec.*0.125);
        dZ = [dZ Real];  % 9,9-Bits (m)
        Bin2Dec = [bin2dec(Temp((57:66)+(i-1)*106)) bin2dec(Temp((108:117)+(i-1)*106))];
        Real = (Bin2Dec > 511).*(Bin2Dec.*(2^-31) - (2^-31)*2^10) + (Bin2Dec <= 511).*(Bin2Dec.*2^-31);
        d_af0 = [d_af0 Real];    % 10,10-Bits (second)
        IODP = [IODP bin2dec(Temp((118:119)+(i-1)*106))];    % 2-Bits
        % Spare = [Spare bin2dec(Temp((120)+(i-1)*106))];    % 1-Bits
    end
end
Output.PRN_Mask = PRN_Mask;
Output.IOD = IOD;
Output.dX = dX;
Output.dY = dY;
Output.dZ = dZ;
Output.dXr = dXr;
Output.dYr = dYr;
Output.dZr = dZr;
Output.d_af0 = d_af0;
Output.d_af1 = d_af1;
Output.t_o = t_o;
Output.IODP = IODP;
clear Temp
