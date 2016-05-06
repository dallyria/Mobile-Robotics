Kgain=-1500;
Dgain=-50;
exitflag=1;
SimTime=.5;

st=111111


for i=1:10
    
ThetaB=StDeg(st);

ac=AcCh(Qr,st);
ThetaE=StDeg(ac);

JointInit=InPosHex(ThetaB);

sim('HexWaRo_StAc_Workspace')
st=ac
i
end
