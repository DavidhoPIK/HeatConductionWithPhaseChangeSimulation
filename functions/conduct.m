function [k] = conduct(T,k_fro,k_nor)
% returns conductivity for given enthalpy and pyhsical data (possibly in matrix form)

liq=(T>0);  % completly liquid

k=liq.*k_nor+(1-liq).*k_fro;
end