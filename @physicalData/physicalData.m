classdef physicalData
    properties
        c_fro  % volumetric heat capacity in frozen state, kJ/m^3/K, x = space, y= time, given each node
        c_nor  % volumetric heat capacity in unfrozend normal state, kJ/m^3/K, given each node
        k_fro  % thermal conductivity in frozen state W/m/K, given each element
        k_nor  % thermal conductivity in frozen state W/m/K, given each element
        L      % volumetric latent heat of fusion, kJ/m^3, given each node
    end
    methods
        function obj= physicalData(c_fro,c_nor,k_fro,k_nor,L)
            obj.c_fro=c_fro;
            obj.c_nor=c_nor;
            obj.k_fro=k_fro;
            obj.k_nor=k_nor;
            obj.L=L
        end

    end
end