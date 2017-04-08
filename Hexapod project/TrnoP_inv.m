% THIS FUNCTION DOES THE TRANSFORMATION FROM HEXAPOD MAIN
% BODY FRAM TO GROUND FRAME 

function Trn_OP_inv=TrnoP_inv(yaw,roll,pitch,OP)
%the robot is moving in Y direction
%yaw: rotationg around Z
%pitch: rotation around x axis
%roll: rotation around y axis
a=-yaw;
b=-roll;
v=-pitch;

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
RxyzT=Rxyz';

Trn_OP_inv=[RxyzT,-RxyzT*OP];
Trn_OP_inv(4,1:4)=[0,0,0,1];
