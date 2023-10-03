%analytical temperature distribution
function T = Tana(t,x,ph,T0,Tini,lambd)
    Tm=0;
    XT = Xana(t,ph,lambd);
    k=ph.k_fro(1,1)/(ph.c_fro(1,1));
    if x<XT
        T= (Tm - T0)/(erf(lambd))*erf(x/(2*(k*t)^(1/2)))+T0;
    elseif x==XT
        T= Tm;
    elseif x>XT
        T=Tini-(Tini-Tm)/erfc(lambd)*erfc(x/(2*(k*t)^(1/2)));
    else
        T=999
    end
end

