function Action=AcCh(Qr,st)

sti=st-111111;
sti=int2str(sti);
Sti=base2dec(sti,3)+1;    
[c,i]=max(Qr(Sti,:));
Aci=dec2base(i-1,3);
Aci=str2num(Aci);
Action=Aci+111111;
