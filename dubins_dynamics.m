function [x,f,g,fh,gh] = dubins_dynamics(paras)

% [x,f,g,fh,gh] = dubins_dynamics(paras) 
% returns the state, f, and g symbolically
% and function handles fh and gh
% 
% State vector X and input vector U is defined as
%   X = [p_x;p_y;theta] = [x1;x2;x3]
%   U = omega = dtheta
%   velocity is assumed to be constant in our case
% Edited by Yi-Hsuan Chen 08/22/2023


%============================== CONSTANT ==================================
v       = paras.v;          % velocity 

%========================= DYNAMICS (ctrlAffine) =========================   
syms p_x p_y theta
x           = [p_x;p_y;theta];
symbolic_f  = [-v*sin(theta); v*cos(theta); 0];
symbolic_g  = [0; 0; 1];

%========================= function handles =========================   
if ~isa(symbolic_f, 'sym')
    f   = sym(symbolic_f);
else
    f   = symbolic_f;
end
if ~isa(symbolic_g, 'sym')
    g   = sym(symbolic_g);
else
    g  = symbolic_g;
end

fh      = matlabFunction(f, 'vars', {x});
gh      = matlabFunction(g, 'vars', {x});
end