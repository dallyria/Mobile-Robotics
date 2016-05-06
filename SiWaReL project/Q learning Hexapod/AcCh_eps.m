function Action=AcCh_eps(Qr,st,eps)

PA=PosAc(st);
a=Qr(Qm(st*ones(size(PA)),PA));
[B,ix]=sort(a);

ixchs=ix(end-floor(length(ix)*eps):end);
ACs=PA(ixchs);
Action=ACs(randi(length(ACs)));