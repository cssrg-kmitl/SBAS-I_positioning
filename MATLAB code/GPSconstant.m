% =========== Constant ================================
GPScons.f1  = 1575.42*10^6;                     %   f1 = 1575.42 MHz (L1)
GPScons.f2  = 1227.60*10^6;                     %   f2 = 1227.60 MHz (L2)
GPScons.we = 7.2921151467e-5;                   %   Earth rotation rate (rad/sec)
GPScons.GM = 3.9860050*10^14;                   %   Earth's universal gravitational parameter (m^3/s^2)
GPScons.c  = 299792458;                         %   Light speed = 299792458 m/s
GPScons.omega_e = 7.2921151467e-5;              %   Earth's rotation rate(rad/sec)
GPScons.lambda1 = 299792458/(1575.42*10^6);     %   Wave length of f1 c/f1
GPScons.lambda2 = 299792458/(1227.60*10^6);     %   Wave length of f2 c/f2
GPScons.k = 9.5196;                             %   Coefficient
GPScons.A = 40.3e16;                            %   Constant value
GPScons.Re = 6371230;                           %   The mean radius of the Earth
GPScons.Rex = 6378136.3;                        %   The radius of the Earth (6378136.3 m)
GPScons.h  = 350*10^3;                          %   Ionospheric height for mapping function (350 km.)
GPScons.ele_cut = 5;                            %   Elevation cut-off