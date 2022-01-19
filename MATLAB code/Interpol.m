function [Output] = Interpol(Wp,Wmax,Wmin,Vmax,Vmin)
%% ======== Interpolation ========================
% Example: [Output] = Interpol(Wp,Wmax,Wmin,Vmax,Vmin)
% Wp = weighting point
% Wmax = Max weight
% Wmin = Min weight
% Vmax = Value of max weight
% Vmin = Value of min weight
% CSSRG-LAB of KMITL, Thailand.
% Version 1 by Somkit Sophan (March 2019)
%% ===============================================
% ===== RTCA DO-229D (A.4.2.4)=====
Output = Vmin + (Vmax-Vmin)*((Wp-Wmin)./(Wmax-Wmin));
