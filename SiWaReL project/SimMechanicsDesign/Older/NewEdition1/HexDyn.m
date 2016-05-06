clc
clear all
OP=[0;0;8];
JointInit=OP(3,1);
yaw=0;
roll=0;
pitch=0;

OP=[-4;0;8];
Q1=[7.5000;12.9904;0];
phi=[0,pi/3,2*pi/3,pi,4*pi/3,5*pi/3];
for i=1:6
Q2(:,i)=(rot(phi(i))*Q1);
end
Q2=Q2([2,1,3],:);
Ltip=Q2;

%% INVERSE KINEMATIC
HexKin

%% DYNAMIC SIMULATOR
sim('HexWaRo_StAc_WS_MD')
BodPOS=BodPos(:,1:3);
% figure(2)
% plot(tout,BodPos(:,3))

figure(3)

n=100;
DyLtip=[LeTi1(n,:);LeTi2(n,:);LeTi3(n,:);LeTi4(n,:);LeTi5(n,:);LeTi6(n,:)];
plot3([DyLtip(:,1);DyLtip(1,1)],[DyLtip(:,2);DyLtip(1,2)],[DyLtip(:,3);DyLtip(1,3)])
hold on
plot3(BodPOS(:,1),BodPOS(:,2),BodPOS(:,3),'or')
axis equal

LeOfSe=mean(DyLtip);
DyLtipNew=DyLtip-ones(6,1)*LeOfSe;
plot3([DyLtipNew(:,1);DyLtipNew(1,1)],[DyLtipNew(:,2);DyLtipNew(1,2)],[DyLtipNew(:,3);DyLtipNew(1,3)],'k')
BodPOSNew=BodPOS-ones(size(BodPOS,1),1)*LeOfSe;
plot3(BodPOSNew(:,1),BodPOSNew(:,2),BodPOSNew(:,3),'sg')
axis equal
hold off
