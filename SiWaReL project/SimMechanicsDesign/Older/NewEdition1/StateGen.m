load ThetaStaCom

St=222222
S6=mod(St,10)
S5=mod(.1*(St-S6),10)
S4=mod(.01*(St-S6-S5*10),10)
S3=mod(.001*(St-S6-S5*10-S4*100),10)
S2=mod(.0001*(St-S6-S5*10-S4*100-S3*1000),10)
S1=mod(.00001*(St-S6-S5*10-S4*100-S3*1000-S2*10000),10)

Theta.t1=ThetaSta.t1(:,floor(S1));
Theta.t2=ThetaSta.t2(:,floor(S2));
Theta.t3=ThetaSta.t3(:,floor(S3));
Theta.t4=ThetaSta.t4(:,floor(S4));
Theta.t5=ThetaSta.t5(:,floor(S5));
Theta.t6=ThetaSta.t6(:,floor(S6));

Kgain=-1500;
Dgain=-50;
JointInit=InPosHex(Theta);
exitflag=1;
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

