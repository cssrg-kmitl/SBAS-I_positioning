function [New_satpos] = Rotation_Z(Sat_pos,dT)
% ===========================================
% Example: [New_satpos] = Rotation_Z(Sat_pos,dT)
% ===========================================
We = 7.2921151467*10^-5;       % Earth rotation rate(rad/sec)
Theta = We * dT;
Rz = [ cos(Theta)    sin(Theta)   0;
      -sin(Theta)    cos(Theta)   0;
       0                0         1];
New_satpos = (Rz*Sat_pos(:))';
end