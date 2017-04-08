% This Program is about to find the Trajectory of Hexapod Walking Robot
% Advanced Robotic - Supervisor: Dr Osguie
% Shahriari Summer of 2012


clc
clear all


Vy=1; %cm/sec
Vx=1; %cm/sec
delta_t=2.5;% sec for each step
steps=20;
Nstep=2;


%   ** This is the part to modify to rectangular hexapod
%% INITIALIZING THE PROBLEM ==============================
l0=2.5;%cm
l1=7.4;%cm
l2=11.4;%cm
lc=l0;
lf=l1;
lt=l2;
global l0 l1 l2 lt lf lc 
P=zeros(3,1);
% h=7.87; % Radios of platform
% P1=P+[h*sin(pi/6);h*cos(pi/6);0];
% P2=P+[h;0;0];
% P3=P+[h*sin(pi/6);-h*cos(pi/6);0];
% P4=P+[-h*sin(pi/6);-h*cos(pi/6);0];
% P5=P+[-h;0;0];
% P6=P+[-h*sin(pi/6);h*cos(pi/6);0];
L = 35
w = 5
P1=P+[-w;L/2;0];
P2=P+[-w;0;0];
P3=P+[-w;-L/2;0];
P4=P+[w;L/2;0];
P5=P+[w;0;0];
P6=P+[w;-L/2;0];
P_MB=[P4,P5,P6,P3,P2,P1];% The edges of the body where limbs begin in main bodie's frame/

% Giving a sample legs tips position:
Q1=[7.5000;12.9904;0];
phi=[0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
for i=1:6
Q2(:,i)=(rot(phi(i))*Q1);

end
%% Setup:

LVy=2*Vy;
LVx=2*Vx;

TotStp=Nstep*steps;
k=1;
Ltip=Q2;
LQ=Q2;
alfa=.5;

% * Sequencing 
if   mod(Nstep,2)==0
    LQ(1:3,1)=Q2(1:3,1)-alfa*[LVx*delta_t;LVy*delta_t;0];
    LQ(1:3,3)=Q2(1:3,3)-alfa*[LVx*delta_t;LVy*delta_t;0];
    LQ(1:3,5)=Q2(1:3,5)-alfa*[LVx*delta_t;LVy*delta_t;0];
elseif mod(Nstep,2)==1
    LQ(1:3,2)=Q2(1:3,2)-alfa*[LVx*delta_t;LVy*delta_t;0];
    LQ(1:3,4)=Q2(1:3,4)-alfa*[LVx*delta_t;LVy*delta_t;0];
    LQ(1:3,6)=Q2(1:3,6)-alfa*[LVx*delta_t;LVy*delta_t;0];
end
    Ltip=LQ;
while Nstep>=1
%% Walking Gait =========================================================
% the position of the frame of main body from ground frame.
  
q=0;
LQ=Ltip;
qStep=TotStp-Nstep*steps;
while q<steps
    qq=qStep+q;
    time=delta_t*qq/steps;
    
    OPPy=Vy*time;
    OPPx=Vx*time;
    OPPz=10;
    LVy=2*Vy;
    LVx=2*Vx;
    w=pi/delta_t;
    StpTime=delta_t*q/steps;
    time_array(k)=time;
    
%   **  This the part to modify the trajectory for swing phase.
    
    DeltaLtipX=LVx*delta_t*.5*(1-cos(w*StpTime));
    DeltaLtipY=LVy*delta_t*.5*(1-cos(w*StpTime));
    DeltaLtipZ=1.5*(1-cos(2*w*StpTime));
    %---------------------------------------------------------

    
    DeltaLtip=[DeltaLtipX;DeltaLtipY;DeltaLtipZ];
    
if mod(Nstep,2)==0
    
    Ltip(1:3,1)=LQ(1:3,1)+DeltaLtip;
    Ltip(1:3,3)=LQ(1:3,3)+DeltaLtip;
    Ltip(1:3,5)=LQ(1:3,5)+DeltaLtip;
    Ltip(1:3,2)=LQ(1:3,2);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
    Ltip(1:3,4)=LQ(1:3,4);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
    Ltip(1:3,6)=LQ(1:3,6);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
elseif mod(Nstep,2)==1

    Ltip(1:3,1)=LQ(1:3,1);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
    Ltip(1:3,3)=LQ(1:3,3);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
    Ltip(1:3,5)=LQ(1:3,5);%+[2*Vx*time;2*Vy*time;3*sin(w*time)];
    Ltip(1:3,2)=LQ(1:3,2)+DeltaLtip;
    Ltip(1:3,4)=LQ(1:3,4)+DeltaLtip;
    Ltip(1:3,6)=LQ(1:3,6)+DeltaLtip;

end
% LQ=Ltip;
    
yaw=0;% Around Z
roll=0;% Around Y
pitch=0;% Around X

%----------------------------------------------------------------------------------------
% ** Modify this part to implement a trajectory for the motion of main body: OPPx OPPy OPPz

OP=[OPPx;OPPy;OPPz];


OPP=OP;
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
q=q+1;
k=k+1;
end
Nstep=Nstep-1;
end


if size(Joints_GF_Array,2)/4==TotStp
    clc
    display('Inverese Kinematic for given trajectory is solved.')
    display('the results of inverse kinematic solution is saved in t_array matrix, every 6 coloumns refers to each step.')
else
    clc 
    display('The robot reached some singular points')
    display('The solved solution in other points are saved into t_array matrix in workspace variables.')
t_array;
end

% Hex_Tra_movie
Hex_gait_movie
HexTra_Plot_Degrees
