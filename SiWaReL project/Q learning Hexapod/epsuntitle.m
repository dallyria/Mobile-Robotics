PosAc(st)
a=Qr(Qm(st*ones(size(PosAc(st))),PosAc(st)))
eps=.00
g=max(a)-(max(a)-min(a))*eps
find(a>g)

[C,ia,ic]= unique(a);
[B,ix]=sort(a);

ixchs=ix(end-floor(length(ix)*eps):end)
a(ixchs)