# SBAS-I_positioning
Publisher: CSSRG-LAB of the King Mongkut's Institute of Technology Ladkrabang (KMITL), Thailand.

Requirements: the MATLAB version of R2019a (64-bits) or R2020a (64-bits).
	       the SBAS system of GAGAN and MSAS based on the SBAS-I standard.

Simulation steps:
       1) Copy a folder of "MATLAB code" to your computer at the path of "C:\".
       2) Open the source code file of "SimulatePerformanceSBAS_I.m" with the MATLAB program.
       3) In the "SimulatePerformanceSBAS_I.m" file, in line 13 - 15 choose the RINEX file types such as the observation file,
          navigation file, and SBAS file of the RINEX version 2 format. Those files must be within the source code folder of "C:\MATLAB code)".
       4) Press button "Run".

Simulation results:
	1. Receiver positioning estimations are consisting of three results.
	   1) The single positioning standard is computed by using the ionospheric delays of the Klobuchar model, which represents "SPS".
	   2) The single positioning estimation is computed by using the SBAS-I corrections of the fast and long-term corrections and ionospheric delays 
	   of the Klobuchar model, which represents "SBAS-I(FL)".
	   3) The single positioning estimation is computed by using the SBAS-I corrections of the fast and long-term corrections and ionospheric delays
	   estimated from the ionospheric grid points (IGP), which represent as "SBAS-I(FLI)".

	2. The "sigma major" and "sigma U" are utilized to calculate the horizontal and vertical protection levels, respectively.
