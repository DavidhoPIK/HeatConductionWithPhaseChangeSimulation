function out = FF(H,T0,ph,h)
% Return the value of the function F from finite element scheme
% ph is a spatial physical vector for a given time with length k
% H is length k+1
z=( (H<0).*H )./ph.c_fro +( (H>ph.L).*(H-ph.L) )./ph.c_nor;
T= [T0; z];
liq=(T>0);
liq1=liq(1:end-1);
liq2=liq(2:end);

k_left  = liq1.*ph.k_nor+(1-liq1).*ph.k_fro;
k_right = liq2.*ph.k_nor+(1-liq2).*ph.k_fro;

G=(k_right.*T(2:end)-k_left.*T(1:end-1))./h;

out=[G;0]-[0;G];
out=out(2:end);


end