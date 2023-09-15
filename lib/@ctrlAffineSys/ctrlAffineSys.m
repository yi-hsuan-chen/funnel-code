% Adopted from Jason Choi https://github.com/HybridRobotics/CBF-CLF-Helper
% Edited by Yi-Hsuan Chen (yhchen91@umd.edu) 09/14/2023

classdef ctrlAffineSys < handle
    properties
        xdim    % state dimension
        udim    % control input dimension
        paras   % system-dependent parameters
        
        % control-affine dynamics: xdot = f(x) + g(x)u
        f       % Drift term as a function handle
        g       % vector field as a function handle
        fpoly   % polynomial dynamics

        % linearized dynamics: xdot = A(t)x + B(t)u
        A       % xdim-by-xdim state matrix as a function handle
        B       % xdim-by-udim input-to-state as a function handle
    end

    methods
        function obj = ctrlAffineSys(paras)
            if nargin < 1
                obj.paras = [];
            else
                obj.paras = paras;
            end
            [x,f,g] = defineSystem(obj,paras);
            obj.init(x,f,g);
        end

        function [x,f,g] = defineSystem(obj,paras)
            % Outputs:  x -- state vector
            %           f -- drift term, expresed symbolically wrt x
            %           g -- control vector field, expresed symbolically wrt x
            x = []; f = []; g = [];
        end

        function dx = dynamics(obj,t,x,u)
            % Inputs:   t -- time, x -- state, u -- control
            % Outputs:  dx -- time derivative of state as a function handle
            dx  = obj.f(x) +  obj.g(x)*u;
        end

        function dx = polynomialdyn(obj,t,x,u)
            % Inputs:   t -- time, x -- state, u -- control
            % Outputs:  dx_poly -- polynomial dynamics as a function handle
            dx  = obj.fpoly(x,u);
        end

        function [A,B] = linearize(obj,x0,u0)
            % Inputs:   x0 -- nominal state traj., u -- nominal control traj.
            % Outputs:  [A,B] -- system matrices as function handles
            A  = @(t) obj.A(x0(t));
            B  = @(t) obj.B(u0(t));
        end

    end
end