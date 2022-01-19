function [Output] = CarrierPhaseSmooth(Obs,Cut_PRN,CodeType,PhaseType,Pn_Pre,L1_Pre,Lambda_Li,alph,t)
%% ================================================
% Objective: To smooth the pseudorange of the code by using the pseudorange of the carrier phase (2.1.4.1.1 of RTCA DO-229D).
% Example: [Output] = CarrierPhaseSmooth(Obs,Cut_PRN,CodeType,PhaseType,Pn_Pr,L1_Pr,Lambda_Li,alph,t)
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (January 2020).
%% ================================================
Output.Pn_Pre = Pn_Pre;
Output.L1_Pre = L1_Pre;
Ep_t = find(Obs.Ep.G == t);
PRN   = Obs.PRN.G(Ep_t);
Ep_t(PRN==Cut_PRN) = [];
PRN(PRN==Cut_PRN) = [];
SOD = (Obs.Date.St(4)*60*60) + (Obs.Date.St(5)*60) + Obs.Date.St(6) + Obs.Ep.G(Ep_t); % Second fo day
Ci = Obs.Data.G(Ep_t,ismember(Obs.Type.G,CodeType))'; % Pseudorange C/A code :L1 (m)
Li = Obs.Data.G(Ep_t,ismember(Obs.Type.G,PhaseType))'; % Carrier phase :L1 (cycle)
Index = isnan(Ci)|isnan(Li); % Matched informations between Ci and Li
if ~isempty(Index)
    Ci(Index) = [];
    Li(Index) = [];
    PRN(Index) = [];
    SOD(Index) = [];
end
% === The carrier phase smoothed method
if sum(Pn_Pre(PRN)==0)>0
    Indx = Pn_Pre(PRN)== 0;
    Pn_Pre(PRN(Indx)) = Ci(Indx);
    L1_Pre(PRN(Indx)) = Li(Indx);
end
Pproj = Pn_Pre(PRN) + (Lambda_Li)*(Li - L1_Pre(PRN));
Pn = alph*Ci + (1-alph)*Pproj; % smoothed pseudorange
Output.Pn = Pn;
Output.PRN = PRN;
Output.Pn_Pre(PRN) = Pn;
Output.L1_Pre(PRN) = Li;
Output.SOD = SOD;
