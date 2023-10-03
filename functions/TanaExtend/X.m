function out= X(t,ph,lam)
    a_l= ph.k_nor(1)/ph.c_nor(1);
    out=2*lam*sqrt(a_l*t);
end