classdef DubinsCar < ctrlAffineSys
    methods
        function [x,f,g] = defineSystem(obj,paras)
            syms p_x p_y theta
            x   = [p_x; p_y; theta];
            f   = [paras.v*sin(theta); paras.v*cos(theta); 0];
            g   = [0;0;1];
        end
    end
end