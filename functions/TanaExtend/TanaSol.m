
function T = TanaSol(t,x,ph,T0,Tini,lambd)
    a_l= ph.k_nor(1)/ph.c_nor(1);
    a_s= ph.k_fro(1)/ph.c_fro(1);
    T=Tini-Tini*erfc(x/(2*sqrt(a_s*t) ))/erfc(lambd*sqrt(a_l/a_s));
end