function T = TanaLiq(t,x,ph,T0,Tini,lambd)
    a_l= ph.k_nor(1)/ph.c_nor(1);
    T=T0-(T0)*erf(x/(2*sqrt(a_l*t)))/erf(lambd);
end