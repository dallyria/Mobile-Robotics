Kgain=-1500;
Dgain=-50;
exitflag=1;
SimTime=.5;

st=121212


for i=1:10
    
ThetaB=StDeg(st);

% ac=AcCh(Qr,st);
ac=AcCh_eps(Qr,st,.05);
ThetaE=StDeg(ac);

JointInit=InPosHex(ThetaB);

sim('EvHexWaRo_StAc_Workspace')
st=ac
i
end
