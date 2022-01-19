R = 1;
Indx_St = 1;
Time_Plot = SODs(:,R)./3600;
%% === Estimated receiver position
figure(1);
cla;
plot(SPS.LLA(Indx_St:end,2),SPS.LLA(Indx_St:end,1)); hold on;
plot(FL.LLA(Indx_St:end,2),FL.LLA(Indx_St:end,1)); hold on;
plot(FLI.LLA(Indx_St:end,2),FLI.LLA(Indx_St:end,1)); hold on;
title('Receiver positioning estimations');
xlabel('Longitude (Degrees)');
ylabel('Latitude (Degrees)');
legend('SPS','SBAS-I(FL)','SBAS-I(FLI)');
grid on;

%% === Standard deviations of the "Mejor" and "U" to compute the HPL and VPL
figure(2);
subplot(211);
cla;
plot(Time_Plot,FLI.SD_Mejor(:,R)); hold on;
title('\sigma_M_e_j_o_r');
xlabel('UTC (h)');
ylabel('Meters');
legend('SBAS-I(FLI)');
xlim([0 24]);
grid on;
subplot(212);
cla;
plot(Time_Plot,FLI.SD_U(:,R)); hold on;
title('\sigma_U');
xlabel('UTC (h)');
ylabel('Meters');
legend('SBAS-I(FLI)');
xlim([0 24]);
grid on;
