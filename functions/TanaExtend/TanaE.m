%analytical temperature distribution
% source: mathematical modeling of melting and freezing processes (V.
% Alexiades) p 48
function T = TanaE(t,x,ph,T0,Tini,lambd)

    if(x<X(t,ph,lambd))
        T= TanaLiq(t,x,ph,T0,Tini,lambd);
    else
        T= TanaSol(t,x,ph,T0,Tini,lambd);
    end
end





