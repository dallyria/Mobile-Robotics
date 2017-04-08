syms yaw pitch roll OPz OPx OPy LTipix LTipiy LTipiz Pix Piy Piz

OP=[OPx;OPy;OPz];
LTipi=[LTipix;LTipiy;LTipiz];
Trn_op=TrnoPe(yaw,pitch,roll,-OP);
X=Trn_op*[LTipi;1]-[Pix;Piy;Piz;1];
X=X(1:3);
simplify(X)

xi=X(1)
yi=X(2)
zi=X(3)





