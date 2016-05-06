% This Program is about to find the Trajectory of Hexapod Walking Robot
% Advanced Robotic - Supervisor: Dr Osguie
% Shahriari Summer of 2012
clc
clear all
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

% Giving a sample legs tips position:
Q1=[7.5000;12.9904;0];
phi=[0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
for i=1:6
Q2(:,i)=(rot(phi(i))*Q1);

end
%% TRAJECTORY:
steps=20;
OPPz=5*ones(1,steps);
OPPx=zeros(1,steps);
OPPy=zeros(1,steps);
YAW=zeros(1,steps);% AROUND Z
ROLL=zeros(1,steps);% AROUND Y
PITCH=zeros(1,steps);% AROUND X
LTIP=[Q2,Q2,Q2,Q2,Q2];
for q=1:steps
    YAW(q)=.28*(q*2*pi/steps);
% %     OPPx(q)=steps*.5*sin(q*2*pi/steps)/5;
%     OPPy(q)=-steps/4*sin(q*2*pi/steps);
%     OPPz(q)=8+4*sin(q*2*pi/steps);
%     OPPz(q)=q+1;
%     ROLL(q)=.3*sin(q*2*pi/steps);
%     PITCH(q)=.3*sin(q*2*pi/steps)

end
OPP=[OPPx;OPPy;OPPz];
k=1;

while k<=steps
%% Configuration =========================================================
OP=OPP(1:3,k); % the position of the frame of main body from ground frame.
Ltip=Q2;
% Ltip(2:3,3)=Ltip(2:3,3)-k*[1,-1]';
yaw=YAW(k);% Around Z
roll=ROLL(k);% Around Y
pitch=PITCH(k);% Around X
%% Define Legs tip in Main body frame====================================

Trn_op=TrnoP(yaw,roll,pitch,-OP);
for i=1:6
    Ltip_MB(1:4,i)=Trn_op*[Ltip(1:3,i);1];
end
Ltip_MB=Ltip_MB(1:3,:);

%% solving inverse Kinematic Problem===========================
% first we should define the X,Y,Z distance of the leg tip from the leg
% frame positioned on the body ( P1,P2...)

Ltip_JF=Ltip_MB-P_MB; % Leg tip in Joint Frame
global Ltip_JF
[t,fval,exitflag]=fsolve('Legs',zeros(3,6));
if exitflag==1
t1=t(1,1:6)';
t2=t(2,1:6)';
t3=t(3,1:6)';
t_array(1:3,6*k-5:6*k)=t;

%% Find All the Robots joint configuration in Main Body Frame=============
% every 3 rows denote the x,y,z coordinates of the joints from body to end
% point on the ground.

for i=1:6
Joints_MB(3*i-2:3*i,1:4)=[P_MB(1,i),...
    P_MB(1,i)+lc*cos(t1(i)),...
    P_MB(1,i)+lc*cos(t1(i))+lf*cos(t2(i))*cos(t1(i)),...
    P_MB(1,i)+cos(t1(i))*(lc+lf*cos(t2(i))+lt*cos(t2(i)+t3(i)));...
    P_MB(2,i),...
    P_MB(2,i)+lc*sin(t1(i)),...
    P_MB(2,i)+sin(t1(i))*(lc+lf*cos(t2(i))),...
    P_MB(2,i)+sin(t1(i))*(lc+lf*cos(t2(i))+lt*cos(t2(i)+t3(i)));...
    P_MB(3,i),...
    P_MB(3,i),...
    P_MB(3,i)-lf*sin(t2(i)),...
    P_MB(3,i)-lf*sin(t2(i))-lt*sin(t3(i)+t2(i))];
end

%% define all the points back into the ground coordinate==================
Trn_po=TrnoP_inv(yaw,roll,pitch,-OP);
for i=1:6
for j=1:4
    J1=Trn_po*[Joints_MB(3*i-2:3*i,j);1];
    J1=J1(1:3);
    Joints_GF(3*i-2:3*i,j)=J1;
end
end
Joints_GF_Array(1:18,4*k-3:4*k)=Joints_GF;% Make an array of all points in each configuration for every 4 coloumns.


end
k=k+1;
end


if size(Joints_GF_Array,2)/4==size(OPP,2)
    clc
    display('Inverese Kinematic for given trajectory is solved.')
    display('the results of inverse kinematic solution is saved in t_array matrix, every 6 coloumns refers to each step.')
else
    clc 
    display('The robot reached some singular points')
    display('The solved solution in other points are saved into t_array matrix in workspace variables.')
t_array;
end