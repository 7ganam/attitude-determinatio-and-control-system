qd=[0;0;0;1];
q0=[0.685;0.695;0.153;0.153];
q0=q0/norm(q0);
w0=[0.53;0.53;0.053]*pi/180;
eul_rad =  quat2eul(w0','XYZ')