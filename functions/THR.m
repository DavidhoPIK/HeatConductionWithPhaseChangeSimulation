function T=THR(x,c_fro,c_nor,L, Fac)
% TH converts enthalphy matrix H to tempreature matrix T based on phyiscal
% data
%
%   F = T(x,ph)
%
%   x is the function input variable (enthalpy)
%   ph is an object of the class physicalData
%
eps=((c_fro+c_nor)/2)./Fac(1);
int= (L.*(1./c_nor)-(L/2)./eps)./(1./c_nor-1./eps);
int_l= 1./eps.*(L./2)./(1./eps-1./c_fro);
    T= (x<int_l).*x./c_fro+ (x>int_l).*(x<int).*((x-L/2)./eps) +(x>int).*(x-L)./c_nor;
    %T= (x<0).*x./c_fro +(x>L).*(x-L)./c_nor;

end