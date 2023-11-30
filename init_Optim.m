function [Vtraj0,rho0,xdots,Phi,intvs,t,x,options] = init_Optim(ts,x0,u0,K,S0,fpoly,k)
% 
% Input:
%   ts      -- time samples
%   x0      -- x0(t) n-by-1 nominal state trajectory.
%   u0      -- u0(t) m-by-1 nominal control input.
%   K       -- K(t) tvlqr control gain matrix.
%   S0      -- S(t) n-by-n matrices (soln to riccati equation).
%   fpoly   -- fpoly(t,x,u) Taylor expanded original system (not closed-loop)
%   k       -- number of time segment
% Output:
%   Vtraj0  --
%   V0      -- initial guess of Lyapunov function in 1-by-(N-1) cells with msspoly
%   rho0    -- intial rho expressed in piecewise polynomial (pp)
%   S0      --
% % %   ui      -- tvlqr CL control input in 1-by-(N-1) cells with msspoly
%   xdots   -- Taylor expanded CL system, corresponding to xbardot in text
%   Phi     -- n-by-n matrix associated with incremental Lyapunov function
%   intvs   -- time intervals in 1-by-(N-1) cell
%   options -- optimization options structure
% Edited by Yi-Hsuan Chen 10/02/2023

%%========================= Optimization settings =========================
options.max_iterations  = 10;
options.converge_tol    = 1e-3;
options.clean_tol       = 1e-8;
options.degL1           = 2;
options.degL0           = 1;
options.degLe           = 2;
options.degV            = 4;
options.xdim            = length(x0(1));
% options.rho0            = 0.1*ones(length(ts),1);
options.gX0             = 0.5*[4 0 0; 0 4 0; 0 0 8];
options.delta_s         = 1e-5;

%%============================ Initialize rho =============================
c           = 2;
rhot        = exp(c*(ts-max(ts))/(max(ts)-min(ts)));
rho0pp      = interp1(ts,rhot,'linear','pp');
% options.gX0 = S0(0)/rhot(1);

%%================ convert piecewise polynomial to msspoly ================
N           = length(ts);
t           = msspoly('t',1);
x           = msspoly('x',options.xdim);        % this is xbar in practice!
u0i         = pp2msspoly(t,spline(ts,u0(ts)));
x0i         = pp2msspoly(t,spline(ts,x0(ts)));
Spp         = interp1(ts,reshape(permute(S0(ts),[3 1 2]),N,options.xdim*options.xdim),'linear','pp');
S0i         = pp2msspoly(t,Spp);
rhoi        = pp2msspoly(t,rho0pp);

%===================== compute (V0,V0dot,rho,rhodot) ======================
% intvs   = cell(1,N-1);
V0      = cell(1,N-1);
rho0    = cell(1,N-1);
S0s     = cell(1,N-1);
xdots   = cell(1,N-1);
Phi     = cell(1,N-1);

for i = 1:N-1     
    intvs{i}    = [ts(i) ts(i+1)]-ts(i);
    S0s{i}      = reshape(S0i(:,i),length(x),length(x));
    % Use matrix 'S' for initial Lyapunov guess
    V0{i}       = x'*S0s{i}*x;
    xdots{i}    = fpoly(t,x+x0i(:,i),u0i(:,i)+K(ts(i))*x)-fpoly(t,x0i(:,i),u0i(:,i));
    V0dot{i}    = diff(V0{i},x)*xdots{i} + diff(V0{i},t);
    Phi{i}      = zeros(length(x),length(x));

    % remove time-dependent
    tt          = linspace(intvs{i}(1),intvs{i}(2),k);
    V0{i}       = subs(V0{i},t,tt);
    V0dot{i}    = subs(V0dot{i},t,tt);
    rho0{i}     = subs(rhoi(:,i),t,tt);
    S0s{i}      = subs(S0s{i},t,tt);
    xdots{i}    = subs(xdots{i},t,tt);
end
Vtraj0.V0       = V0;
Vtraj0.V0dot    = V0dot;
Vtraj0.S0s      = S0s;
end