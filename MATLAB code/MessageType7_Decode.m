function [Output] = MessageType7_Decode(Type7)
%% ==============================================================
% Objective: To convert the message types 7 to the fast correction degradation factor (A.4.4.5 of RTCA DO-229D).
% Example: [Output] = MessageType7_Decode(Type7).
% Type7 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================================
% === Constant ===
c_ai = [0 0.00005 0.00009 0.00012 0.00015 0.00020 0.00030 0.00045 ...
    0.00060 0.00090 0.00150 0.00210 0.00270 0.00330 0.00460 0.00580]; % Fast Corrections Degradation Factor (ai) m/s^2
c_Ifc1 = [180 180 153 135 135 117 99 81 63 45 45 27 27 27 18 18]; % User Time-Out Interval for fast corrections-seconds En Route through LNAV Approach 
c_Ifc2 = [120 120 102 90 90 78 66 54 42 30 30 18 18 18 12 12]; % User Time-Out Interval for fast corrections-seconds LNAV/VNAV, LPV, LP Approach
c_s = [60 60 51 45 45 39 33 27 21 15 15 9 9 9 6 6]; % Maximum Fast Correction Update Interval (second)
% === Parameter variables ===
PRN = 51; % PRN number
ai = nan(1,PRN);    % Fast corrections degradation indicator
ai2 = nan(1,PRN);   % Fast corrections degradation fector (m/s^2)
Ifc1 = nan(1,PRN);  % User time-out interval for fast corrections-seconds En-route through LNAV approach
Ifc2 = nan(1,PRN);  % User time-out interval for fast corrections-seconds LNAV/VNAV,LPV, LP approach
t_Update = nan(1,PRN);  % Maximum fast correction update interval (seconds)
% ===========================
Char_hex = char(Type7)';    % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Output.t_lat = bin2dec(Temp(15:18)); % System latency for 4-Bits
Output.IODP = bin2dec(Temp(19:20));  % IODP status for 2-Bits
% Spare = bin2dec(Temp(21:22)); % Spare for 2-Bits

for i = 1:51 % For 51x4-Bits
    ai(i) = bin2dec(Temp((23:26)+(i-1)*4)); % Degradation factor indicator
    ai2(i) = c_ai(ai(i)+1);   % Degradation factor (m/s^2)
    Ifc1(i) = c_Ifc1(bin2dec(Temp((23:26)+(i-1)*4))+1);  % User Time-Out Interval for NPA mode (s)
    Ifc2(i) = c_Ifc2(bin2dec(Temp((23:26)+(i-1)*4))+1);  % User Time-Out Interval for PA mode (s)
    t_Update(i) = c_s(bin2dec(Temp((23:26)+(i-1)*4))+1); % Maximum Update time(s)
end
Output.ai = ai;
Output.ai2 = ai2;
Output.Ifc1 = Ifc1;
Output.Ifc2 = Ifc2;
Output.t_Update = t_Update;
clear Temp;
