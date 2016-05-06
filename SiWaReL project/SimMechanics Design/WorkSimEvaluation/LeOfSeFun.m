% This function finds the legs those move the robot's body by pushing the
% ground. And the robot's translational movement is approaximated by
% average movement of the legs moved on the ground.

function LegOff=LeOfSeFun(LeTi1,LeTi2,LeTi3,LeTi4,LeTi5,LeTi6)

Be=3; % Begining index, an initial deadline for simulation has been considered.


DyLtipB=[LeTi1(Be,:);LeTi2(Be,:);LeTi3(Be,:);LeTi4(Be,:);LeTi5(Be,:);LeTi6(Be,:)];
DyLtipE=[LeTi1(end,:);LeTi2(end,:);LeTi3(end,:);LeTi4(end,:);LeTi5(end,:);LeTi6(end,:)];

[rb,cb]=find(DyLtipB(:,3)>.015);
[re,ce]=find(DyLtipE(:,3)>.015);

MusDel=[rb;re];
MusDelRow=unique(MusDel);

DyLtipBt=DyLtipB;
DyLtipBt(MusDelRow,:)=[];

DyLtipEt=DyLtipE;
DyLtipEt(MusDelRow,:)=[];

DyLtip=DyLtipEt-DyLtipBt;

LegOff=mean(DyLtip,1);

