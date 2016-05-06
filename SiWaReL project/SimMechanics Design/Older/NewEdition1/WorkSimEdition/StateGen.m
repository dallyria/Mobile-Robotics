clear all
Qr=zeros(729);

st=111111;
ThetaB=StDeg(st);

ac=222111;

ThetaE=StDeg(ac);

% Theta=ThetaB;


Qr(Qm(st,ac))=10;

Kgain=-1500;
Dgain=-50;
JointInit=InPosHex(ThetaB);
exitflag=1;
SimTime=.5;
sim('HexWaRo_StAc_Workspace')


% BodPOS=BodPos(:,1:3);
% figure(2)
% plot(tout,BodPos(:,3))
% figure(3)
% 
% 
% DyLtip=[LeTi1(end,:);LeTi2(end,:);LeTi3(end,:);LeTi4(end,:);LeTi5(end,:);LeTi6(end,:)];
% plot3([DyLtip(:,1);DyLtip(1,1)],[DyLtip(:,2);DyLtip(1,2)],[DyLtip(:,3);DyLtip(1,3)])
% hold on
% plot3(BodPOS(:,1),BodPOS(:,2),BodPOS(:,3),'or')
% axis equal
% 
% LeOfSe=mean(DyLtip);
% DyLtipNew=DyLtip-ones(6,1)*LeOfSe;
% plot3([DyLtipNew(:,1);DyLtipNew(1,1)],[DyLtipNew(:,2);DyLtipNew(1,2)],[DyLtipNew(:,3);DyLtipNew(1,3)],'k')
% BodPOSNew=BodPOS-ones(size(BodPOS,1),1)*LeOfSe;
% plot3(BodPOSNew(:,1),BodPOSNew(:,2),BodPOSNew(:,3),'sg')
% axis equal
% hold off
BodyMov=LeOfSeFun(LeTi1,LeTi2,LeTi3,LeTi4,LeTi5,LeTi6)
BodyTilt=BoOriFun(BodPos(:,4:6))