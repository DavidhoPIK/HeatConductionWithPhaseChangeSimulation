
% function to determine parmeter for analytical solution
function Q = evalLam(x,ph,T0,Tini)
Tm=0;
K_s =ph.k_fro(1,1);
K_l=ph.k_nor(1,1);
C_s=ph.c_fro(1,1);
k_s = K_s/(C_s);
k_l = K_l/(C_s);
f=@(x) exp(-x^2)/erf(x)- K_l/K_s*...
    (sqrt(k_s)*(Tini-Tm)*exp(-k_s*x^2/k_l))/...
    (sqrt(k_l)*(Tm-T0)*erfc(x*sqrt(k_s/k_l)))-...
    (x*ph.L(1,1)*sqrt(pi))/(C_s*(Tm-T0));
Q=f(x);
end