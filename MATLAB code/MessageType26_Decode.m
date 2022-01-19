function [Output] = MessageType26_Decode(Type26)
%% ==============================================================
% Objective: To convert the message types 26 to the ionospheric delay corrections(A.4.4.10 of RTCA DO-229D).
% Example: [Output] = MessageType26_Decode(Type26).
% Type26 is an input as the Hex number in the string type (cell size 1x1).
% Output is a structure data type as the matrix of double values.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019).
%% ==============================================================
GIVEi = [0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3 3.6 4.5 6 15 45 NaN];    % Grid Ionospheric Vertical Error Indicator(GIVEI) Constant
Sig2_GIVEi = [0.0084 0.0333 0.0749 0.1331 0.2079 0.2994 0.4075 0.5322...
    0.6735 0.8315 1.1974 1.8709 3.3260 20.7870 187.0826 NaN];       % Sigma^2 of Grid Ionospheric Vertical Error Indicator Constant
% === Data buff ===
Temp1 = NaN(15,1);
Temp2 = NaN(15,2);
Temp3 = NaN(6,1);
Temp4 = NaN(6,2);
Temp5 = NaN(5,1);
Temp6 = NaN(5,2);
Temp7 = NaN(12,1);
Temp8 = NaN(12,2);
% =================
Char_hex = char(Type26)';   % A cell array of character vectors to a row of the character array
Bin = dec2bin(hex2dec(Char_hex),4)';
Temp = Bin(:)';
Band_numb = bin2dec(Temp(15:18));    % Band Numbers
Block_ID = bin2dec(Temp(19:22));    % Block ID
Output.Band_numb = Band_numb;
Output.Block_ID = Block_ID;
if (Band_numb<8)
    if (Block_ID<13)
        for i = 1:15
            Temp1(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp2(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp1;
        GIVE = Temp2(:,1);
        Sig2_GIVE = Temp2(:,2);
    else
        for i = 1:6
            Temp3(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp4(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp3;
        GIVE = Temp4(:,1);
        Sig2_GIVE = Temp4(:,2);
    end
elseif (Band_numb==8)
    if (Block_ID<13)
        for i = 1:15
            Temp1(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp2(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp1;
        GIVE = Temp2(:,1);
        Sig2_GIVE = Temp2(:,2);
    else
        for i = 1:5
            Temp5(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp6(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp5;
        GIVE = Temp6(:,1);
        Sig2_GIVE = Temp6(:,2);
    end
elseif (Band_numb>8)
    if (Block_ID<12)
        for i = 1:15
            Temp1(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp2(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp1;
        GIVE = Temp2(:,1);
        Sig2_GIVE = Temp2(:,2);
    else
        for i = 1:12
            Temp7(i) = bin2dec(Temp((23:31)+(i-1)*13))*0.125; % IGP Vertical Delay Estimaste (meters)
            Temp8(i,:) = [GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1) Sig2_GIVEi(bin2dec(Temp((32:35)+(i-1)*13))+1)]; % Grid Ionospheric Vertical Error Indicator
        end
        IGP_VDE = Temp7;
        GIVE = Temp8(:,1);
        Sig2_GIVE = Temp8(:,2);
    end
end
IODI = bin2dec(Temp(218:219)); % IODI (Issue of Data Ionospheric) status
% Spare = bin2dec(Temp(220:226));
Output.IGP_VDE = IGP_VDE;
Output.GIVE = GIVE;
Output.Sig2_GIVE = Sig2_GIVE;
Output.IODI = IODI;
clear Temp1 Temp2 Temp3 Temp4 Temp5 Temp6 Temp7 Temp8
