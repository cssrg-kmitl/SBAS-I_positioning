function [Output] = MessageType2to5_Decode(TypeX)
%% ==============================================
% Objective: To convert the message types 2-5 to the fast corrections (A.4.4.3 of RTCA DO-229D).
% Example: [Output] = MessageType2to5_Decode(TypeX).
% TypeX is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================
UDREi = [0.75 1 1.25 1.75 2.25 3 3.75 4.5 5.25 6 7.5 15 50 150 NaN NaN]';    % UDREi constant table (meters)
Sigma2 = [0.052 0.0924 0.1444 0.283 0.4678 0.8315 1.2992 1.8709 2.5465...
    3.3260 5.1968 20.787 230.9661 2078.695 NaN NaN]';   % Sigma2_UDREi constant table (meters^2)
Char_hex = char(TypeX)';    % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Output.IODF = bin2dec(Temp(15:16)); % IOPF status numbers
Output.IODP = bin2dec(Temp(17:18)); % IODP status numbers
for i = 1:13
    Bin2Dec = bin2dec(Temp((19:30)+(i-1)*12));
    if (Bin2Dec <= 2047)
        PRC(i) = Bin2Dec*0.125; % Convert binary bits to Fast corrections (Positive values)
    else
        PRC(i) = Bin2Dec*0.125 - 256*2; % Convert binary bits to Fast corrections (Negative values)
    end
    Index = bin2dec(Temp((175:178)+(i-1)*4))+1;
    UDRE(i) = UDREi(Index); % UDREi values
    Sig2_UDRE(i) = Sigma2(Index);   % Sigma^2 of UDREi values
end
Output.PRC = PRC;
Output.UDRE = UDRE;
Output.Sig2_UDRE = Sig2_UDRE;
clear Temp;
