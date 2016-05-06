function NewAc=SofMaxAc(cu_st)

sap=Qm(cu_st*ones(64,1),PosAc(cu_st));%All possible future actions;
exp(Qr(sap))/