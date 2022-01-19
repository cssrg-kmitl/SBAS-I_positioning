function [Output] = MessageType1_Decode(Type1)
%% ====================================================
% Objective: To convert the message type 1 to the PRN mask assignments (A.4.4.2 of RTCA DO-229D).
% Example: [Output] = MessageType1_Decode(Type1).
% Type1 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ====================================================
Char_hex = char(Type1)';    % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';	% Convert characters to decimal numbers
Temp = Bin(:)';
bPRNmask = double(Temp(15:224))>48;	% 210-PRN Slot of Binary bits
Output.IODP = bin2dec(Temp(225:226));  % IODP status numbers
Output.PRNcode = bPRNmask;	% PRN code numbers
clear Temp;
Temp = find(bPRNmask);
Output.PRN_Mask(1:length(Temp)) = find(Temp); % PRN mask numbers
clear Temp;
