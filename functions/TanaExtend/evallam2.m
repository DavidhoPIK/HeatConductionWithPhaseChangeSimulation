function out = evallam2(x,ph,T0,Tini)
    St_l = ph.c_nor(1)*(T0)/ph.L(1);
    St_s = ph.c_fro(1)*(-Tini)/ph.L(1);
    a_l= ph.k_nor(1)/ph.c_nor(1);
    a_s= ph.k_fro(1)/ph.c_fro(1);
    v=sqrt(a_l/a_s);
    out=(St_l/(exp(x*x) * erf(x)))- (St_s / (v *exp(v^2* x^2) *erfc(v*x))) - x*sqrt(pi);
end