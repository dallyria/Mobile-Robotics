% this function finds the indexes of Q matrix using the state and
% action.

function Y=Qm(st,ac)
sti=st-111111;
sti=int2str(sti);

aci=ac-111111;
aci=int2str(aci);

Sti=base2dec(sti,3)+1;
Aci=base2dec(aci,3)+1;

Y=sub2ind([729,729],Sti,Aci);