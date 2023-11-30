function dx = dubins_model_controlled(x,u)
% Inputs:
%    x - state vector
%    u - input vector (here: reference trajectory)
%
% Outputs:
%    f - time-derivative of the state vector
%
% state x = [X,Y,psi]
% input u = [delta,omega_f,omega_r]

% state
X   = x(1); %#ok<NASGU>
Y   = x(2); %#ok<NASGU>
psi = x(3);

%control action
% inputs are values of the state feedback matrix K, the reference state X0,
% and the reference input U0
K   = [u(1) u(2) u(3)];
X0  = [u(4); u(5); u(6)];
U0  = [u(7)];
ufb = U0+K*(x-X0);

% dynamics
v           = 1;
dx(1,1)     = -v*sin(psi);
dx(2,1)     = v*cos(psi);
dx(3,1)     = ufb;
end