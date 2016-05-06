Kgain=-1500;
Dgain=-50;
exitflag=1;
SimTime=.5;

st=222223;


for i=0:5
    
ThetaB=StDeg(st);

ac=AcCh(Qr,st);
ThetaE=StDeg(ac);

JointInit=InPosHex(ThetaB);

sim('HexWaRo_StAc_Workspace')
st=ac
end
