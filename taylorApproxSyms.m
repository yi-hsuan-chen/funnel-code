function polysyms = taylorApproxSyms(fh,gh,order)
% calculates the polynomial dynamics of a non polynomial system
% by taylor-expanding around the nominal trajectory to a certain degree.
% Input:
%   fh     -- drift term, expressed symbolically wrt x.
%   gh     -- control vector fields, expressed symbolically wrt x.
%   order  -- specify truncation order.
% Output:
%   polysyms  -- the function handle of @(t,x,u) polynomial dynamics.
% Edited by Yi-Hsuan Chen 08/25/2023

syms t u
x           = sym('x',[3 1]);
fcl         = fh(x) + gh(x)*u;
poly        = vpa(taylor(fcl,x,'Order',order));
polysyms    = matlabFunction(poly, 'vars', {t,x,u});
end