clear all

alpha=.5;
lambda=.9;
Qr=zeros(729);
Rew_mat=zeros(729);
Fuzzy_Reward=readfis('Reward');
Kgain=-1500;
Dgain=-50;
exitflag=1;
SimTime=.5;

st=111111;

tic;
for iter=1:10000

ThetaB=StDeg(st);

PossibAc=PosAc(st);
ac=PossibAc(randi([1,64]));

ThetaE=StDeg(ac);

JointInit=InPosHex(ThetaB);

sim('HexWaRo_StAc_Workspace')

BodyMov=LeOfSeFun(LeTi1,LeTi2,LeTi3,LeTi4,LeTi5,LeTi6);
BodyTilt=BoOriFun(BodPos(:,4:6));

Rew=evalfis([100*BodyMov,BodyTilt],Fuzzy_Reward);% Displacement in Fuzzy Reward is assigned in (cm)

sa=Qm(st,ac);
% [st_i,ac_i]=ind2sub([729,729],sa);

cu_st=ac;
FuAc=PosAc(cu_st);
sap=Qm(cu_st*ones(64,1),FuAc);%All possible future actions;



Rew_mat(sa)=Rew;
Qr(sa)=alpha*(Rew+lambda*max(max(Qr(sap))))+(1-alpha)*Qr(sa);

st=ac
end

EvalTime=toc;
% save LearnRes Qr st