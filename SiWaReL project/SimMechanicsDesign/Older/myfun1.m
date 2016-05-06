function [ta1,ta2,ta3,ta4,ta5,ta6,exitflag] = myfun1(OP,roll,pitch,yaw)
%#codegen

% THIS PROGRAM IS GOING TO TAKE ALL THE GROUND COORDINATES TO MOVING
% PLATFORM FRAME AND CALCULATE INVERSE KINEMATIC THERE.
% ADVANCED ROBOTIC PROJECT - Supervisor Dr Osguie
%  - Shahriari Summer 2012

%% INITIALIZING THE PROBLEM ==============================
l0=2.5;%cm
l1=7.4;%cm
l2=11.4;%cm
lc=l0;
lf=l1;
lt=l2;
global l0 l1 l2 lt lf lc 
P=zeros(3,1);
h=7.87; % Radios of platform
P1=P+[h*sin(pi/6);h*cos(pi/6);0];
P2=P+[h;0;0];
P3=P+[h*sin(pi/6);-h*cos(pi/6);0];
P4=P+[-h*sin(pi/6);-h*cos(pi/6);0];
P5=P+[-h;0;0];
P6=P+[-h*sin(pi/6);h*cos(pi/6);0];
P_MB=[P1,P2,P3,P4,P5,P6];% The edges of the body where limbs begin in main bodie's frame/
P_MB=P_MB([2,1,3],:);
% Giving a sample legs tips position:
Q1=[7.5000;12.9904;0];
phi=[0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
for i=1:6
RT=[cos(phi(i)),sin(phi(i)),0;-sin(phi(i)),cos(phi(i)),0;0,0,1];
Q2(:,i)=(RT*Q1);end
Q2=Q2([2,1,3],:);
%% Configuration =========================================================
% OP=[0;0;6]; % the position of the frame of main body from ground frame.
Ltip=Q2;
% yaw=0;% Around Z
% roll=0;% Around X
% pitch=0;% Around Y
%% Define Legs tip in Main body frame====================================
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


Trn_op=Trn_OP;
for i=1:6
    Ltip_MB(:,i)=Trn_op*[Ltip(:,i);1];
end
Ltip_MB=Ltip_MB(1:3,:);

%% solving inverse Kinematic Problem===========================
% first we should define the X,Y,Z distance of the leg tip from the leg
% frame positioned on the body ( P1,P2...)

Ltip_JF=Ltip_MB-P_MB; % Leg tip in Joint Frame
% global Ltip_JF
% [t,fval,exitflag]=fsolve('Legs',zeros(3,6));
for i=1:6    
    x=Ltip_JF(1,i);
    y=Ltip_JF(2,i);
    z=Ltip_JF(3,i);
       
        t1(i)=atan2(y,x);
        c=((x-l0*cos(t1(i)))^2+(y-l0*sin(t1(i)))^2+z^2)^.5;
        B=acos((-l2^2+l1^2+c^2)/(2*l1*c));
        t2(i)=(asin(-z/c)-B);
        C1=pi/2-B;
        h=l1*sin(B);
        C2=acos(h/l2);
        t3(i)=(pi-C1-C2);
          
end
ExFl=imag(t1)+imag(t2)+imag(t3);
ExFl=mean(ExFl);
if ExFl==0
exitflag=1;
else
    exitflag=-1;
end
t1=t1*180/pi;
t2=t2*180/pi;
t3=t3*180/pi;
ta1=[t1(1),t2(1),t3(1)]';
ta2=[t1(2),t2(2),t3(2)]';
ta3=[t1(3),t2(3),t3(3)]';
ta4=[t1(4),t2(4),t3(4)]';
ta5=[t1(5),t2(5),t3(5)]';
ta6=[t1(6),t2(6),t3(6)]';


