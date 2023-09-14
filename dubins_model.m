function XDOT = dubins_model(X,U,paras)

% [XDOT] = dubins_model(X,U) returns the state derivatives, XDOT,
% given the current state and inputs, X and U.
%
% State vector X and input vector U is defined as
%   X = [p_x;p_y;theta] = [x1;x2;x3]
%   U = omega = dtheta
%   velocity is assumed to be constant in our case
% Edited by Yi-Hsuan Chen 08/22/2023

%====================== STATE AND CONTROL VECTOR ==========================

%%... state variables
p_x     = X(1,:);           % x position in globle frame
p_y     = X(2,:);           % y position in globle frame
theta   = X(3,:);           % yaw angle (orientation)
%%... control inputs
u       = U(1,:);           % steering rate

%============================== CONSTANT ==================================
v       = paras.v;          % velocity 

%========================= STATE DERIVATIVES =========================   
XDOT    = [v*sin(theta); v*cos(theta); u];
end