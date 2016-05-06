% THIS FUNCTION DOES THE TRANSFORMATION FROM GROUND FRAME TO HEXAPOD MAIN
% BODY FRAM

function Trn_OP=TrnoPe(yaw,pitch,roll,OP)
%the robot is moving in Y direction
%yaw: rotationg around Z
%pitch: rotation around x axis
%roll: rotation around y axis
a=-yaw;
b=-pitch;
v=-roll;

Rz=[cos(a),-sin(a),0;...
    sin(a),cos(a),0;...
    0,0,1];
Ry=[cos(b),0,sin(b);...
    0,1,0;...
    -sin(b),0,cos(b)];
Rx=[1,0,0;...
    0,cos(v),-sin(v);
    0,sin(v),cos(v)];
Rxyz=Rz*Ry*Rx;

Trn_OP=[Rxyz,OP];
Trn_OP(4,1:4)=[0,0,0,1];
