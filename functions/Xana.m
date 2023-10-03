
% analytical position of freezing front
function X = Xana(t,ph,lambd)
k_s = ph.k_fro(1,1)/(ph.c_fro(1,1));
X= 2*lambd*sqrt(k_s * t);
end
