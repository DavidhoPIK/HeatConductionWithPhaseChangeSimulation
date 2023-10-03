function T=TH(x,c_fro,c_nor,L)
% TH converts enthalphy matrix H to tempreature matrix T based on phyiscal
% data
%
%   F = T(x,ph)
%
%   x is the function input variable (enthalpy)
%   ph is an object of the class physicalData
%
    T= (x<0).*x./c_fro +(x>L).*(x-L)./c_nor;

end