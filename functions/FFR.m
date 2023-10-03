function out = FFR(H,ph,h, rfac)
% Return the value of the function F from finite element scheme
% ph is a spatial physical vector for a given time with length k
% H is length k+1

k_left = conduct(H(1:end-1),ph.k_fro,ph.k_nor,ph.L(1:end-1));
T_left = THR(H(1:end-1),ph.c_fro(1:end-1),ph.c_nor(1:end-1),ph.L(1:end-1),rfac );

k_right = conduct(H(2:end),ph.k_fro,ph.k_nor,ph.L(2:end));
T_right = THR(H(2:end),ph.c_fro(2:end),ph.c_nor(2:end),ph.L(2:end),rfac );

G=(k_right.*T_right-k_left.*T_left)./h;

out=[G;0]-[0;G];


end