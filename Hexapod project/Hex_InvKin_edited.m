% THIS PROGRAM IS GOING TO TAKE ALL THE GROUND COORDINATES TO MOVING
% PLATFORM FRAME AND CALCULATE INVERSE KINEMATIC THERE.
% ADVANCED ROBOTIC PROJECT - Supervisor Dr Osguie
%  - Shahriari Summer 2012
tic

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
P_MB=P_MB([2,1,3],:);
% Giving a sample legs tips position:
Q1=[7.5000;12.9904;0];
phi=[0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
for i=1:6
Q2(:,i)=(rot(phi(i))*Q1);
end
Q2=Q2([2,1,3],:);
%% Configuration =========================================================
OP=[0;0;6]; % the position of the frame of main body from ground frame.
Ltip=Q2;
yaw=0;% Around Z
roll=0;% Around X
pitch=0;% Around Y
%% Define Legs tip in Main body frame====================================

Trn_op=TrnoPe(yaw,pitch,roll,-OP);
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

if exitflag==1
% t1=t(1,1:6)';
% t2=t(2,1:6)';
% t3=t(3,1:6)';


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
Trn_po=TrnoPe_inv(yaw,pitch,roll,-OP);
for i=1:6
for j=1:4
    J1=Trn_po*[Joints_MB(3*i-2:3*i,j);1];
    J1=J1(1:3);
    Joints_GF(3*i-2:3*i,j)=J1;
end
end

BJ=[[Joints_GF(1,1),Joints_GF(4,1),Joints_GF(7,1),...
    Joints_GF(10,1),Joints_GF(13,1),Joints_GF(16,1),Joints_GF(1,1)];...
    [Joints_GF(2,1),Joints_GF(5,1),Joints_GF(8,1),...
    Joints_GF(11,1),Joints_GF(14,1),Joints_GF(17,1),Joints_GF(2,1)];...
    [Joints_GF(3,1),Joints_GF(6,1),Joints_GF(9,1),...
    Joints_GF(12,1),Joints_GF(15,1),Joints_GF(18,1),Joints_GF(3,1)]];% Body joints in ground frame
%% Plot Hexapod===========================================================
figure(2)
for i=1:6
    plot3(Joints_GF(3*i-2,:),Joints_GF(3*i-1,:),Joints_GF(3*i,:),'-bo','Linewidth',2,...
        'MarkerEdgeColor','r',...
                'MarkerFaceColor','c',...
                'MarkerSize',6,'linewidth',7)
    hold on
end

for i=1:6
fill3([Joints_GF(1,1),Joints_GF(4,1),Joints_GF(7,1),...
    Joints_GF(10,1),Joints_GF(13,1),Joints_GF(16,1),Joints_GF(1,1)],...
    [Joints_GF(2,1),Joints_GF(5,1),Joints_GF(8,1),...
    Joints_GF(11,1),Joints_GF(14,1),Joints_GF(17,1),Joints_GF(2,1)],...
    [Joints_GF(3,1),Joints_GF(6,1),Joints_GF(9,1),...
    Joints_GF(12,1),Joints_GF(15,1),Joints_GF(18,1),Joints_GF(3,1)],'c','Linewidth',3)
    
end
hold off
axis equal
Grid on
text(BJ(1,1),BJ(2,1),BJ(3,1)+1,'#1')
text(BJ(1,2),BJ(2,2),BJ(3,2)+1,'#2')
text(BJ(1,3),BJ(2,3),BJ(3,3)+1,'#3')
text(BJ(1,4),BJ(2,4),BJ(3,4)+1,'#4')
text(BJ(1,5),BJ(2,5),BJ(3,5)+1,'#5')
text(BJ(1,6),BJ(2,6),BJ(3,6)+1,'#6')

xlabel ' X(cm) '%changed should be corrected completely
ylabel ' Y(cm) '
zlabel ' Z(cm) '   
title 'Hexapod Walking Robot'
clc
display('Inverse Kinematic Problem is solved; Leg joints angles found.')

display('-------------------------------------------------------------')
disp('      t1       t2         t3')
disp([t1,t2,t3]*180./pi)
display('*The angles are defined in the Main body frame coordinates.')

else
    clc
    display('Singularity happened. No solution found.')
end
toc