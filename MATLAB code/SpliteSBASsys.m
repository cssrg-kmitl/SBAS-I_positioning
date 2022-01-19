function [Output] = SpliteSBASsys(SBAS_File)
%% ==========================================================
% Objective: To splite the each SBAS system and SBAS parameter to a structure data type.
% Exsample: [Output] = SpliteSBASsys(SBAS_File).
% SBAS_File is an input file as the string type obtained from the SBAS systemes.
% Output is a structure data type as the matrix double and cell string types.
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (18 March 2019).
%% ==========================================================
Fid = fopen(SBAS_File); 
C = textscan(Fid,'%s','Delimiter','\n'); % read file (string)
fclose(Fid); % close file
if ~isempty(C)
    Data_split = regexp(C{1,1},'\s*','split'); % Split the column of cell string types
    clear C; % clare variable
    Cell_str = vertcat(Data_split{:}); % Concatenate arrays vertically
    clear Data_split; % clare variable
    PRNs = str2double(Cell_str(:,3)); % Convert strings to double precision values
    N_PRNs = unique(PRNs); % Unique values in PRNs array
    PRN_GAGAN = [127;128;132]; % PRN numbers for GAGAN SBAS systems
    PRN_MSAS = [129;137]; % PRN numbers for MSAS SBAS systems
    
    GAGAN_PRN = []; % Creat empty variable
    for i = 1:length(PRN_GAGAN) % Check the PRN number of the GAGAN SBAS system
        if sum(N_PRNs==PRN_GAGAN(i))
            GAGAN_PRN = [GAGAN_PRN; N_PRNs(N_PRNs==PRN_GAGAN(i))];
        end
    end
    Output.GAGAN_PRN = GAGAN_PRN; % Return output of the GAGAN_PRN variable
    
    MSAS_PRN = []; % Creat empty variable
    for i = 1:length(PRN_MSAS) % Check the PRN number of the MSAS SBAS system
        if sum(N_PRNs==PRN_MSAS(i))
            MSAS_PRN = [MSAS_PRN; N_PRNs(N_PRNs==PRN_MSAS(i))];
        end
    end
    Output.MSAS_PRN = MSAS_PRN; % Return output of the MSAS_PRN variable
    
    if ~isempty(GAGAN_PRN)
        Output.GAGAN_Data = cell(length(GAGAN_PRN),1);
        % === Splite the GAGAN SBAS parameters for the each PRN number
        for n = 1:length(GAGAN_PRN)
            InxType = find(PRNs == GAGAN_PRN(n));
            Output.GAGAN_Data{n} = [Cell_str(InxType,1:4) Cell_str(InxType,6)];
        end
    else
        Output.GAGAN_Data = [];
    end
    
    if ~isempty(MSAS_PRN)
        Output.MSAS_Data = cell(length(MSAS_PRN),1);
        % === Splite the MSAS SBAS parameters for the each PRN number 
        for n = 1:length(MSAS_PRN)
            InxType = find(PRNs == MSAS_PRN(n));
            Output.MSAS_Data{n} = [Cell_str(InxType,1:4) Cell_str(InxType,6)];
        end
    else
        Output.MSAS_Data = [];
    end

else
    error(['The file of' SBAS_File 'is empty']);
end

clear Cell_str;
