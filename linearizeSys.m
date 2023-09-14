function [A,B]  = linearizeSys(x,f,g,x0,u0)
% Input:
%   x     -- symbolic state vector
%   f     -- drift term, expressed symbolically wrt x.
%   g     -- control vector fields, expressed symbolically wrt x.
%   x0    -- x0(t) function n-by-1 double to linearize around.
%   u0    -- u0(t) function m-by-1 double to linearize around.
% Output:
%   A     -- A(t) n-by-n Jacobian of xdot w.r.t. x at x0(t).
%   B     -- B(t) n-by-m Jacobian of xdot w.r.t. u at u0(t).

syms u
xdot    = f + g*u;
Atemp   = jacobian(xdot,x);
Btemp   = g;
A       = @(t) double(subs(Atemp,x,x0(t)));
B       = @(t) double(subs(Btemp,u,u0(t)));
end