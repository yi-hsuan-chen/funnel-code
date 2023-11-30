function [Vtraj0,rho0pp,ui,intvs,t,x,options] = InitializeOptim(ts,x0,u0,K,S0)
% 
% Input:
%   ts      -- time samples
%   fpoly   -- fpoly(t,x,u) Taylor expanded original system (not closed-loop)
%   x0      -- x0(t) n-by-1 nominal state trajectory.
%   u0      -- u0(t) m-by-1 nominal control input.
%   K       -- K(t) tvlqr control gain matrix.
%   S0      -- S(t) n-by-n matrices (soln to riccati equation).
% Output:
%   Vtraj0  -- structure including V0, V0k(dV0/dt), S0s in 1-by-(N-1) cells with msspoly
%   rho0pp  -- intial rho expressed in piecewise polynomial (pp)
% % %   rho0    -- struture including rho and rhodot in 1-by-(N-1) cells with doubles
%   Phi     -- initial Phi in 1-by-(N-1) cell with n-by-n matrices 
%   ui      -- tvlqr CL control input in 1-by-(N-1) cells with msspoly
%   intvs   -- time intervals in 1-by-(N-1) cell
%   options -- optimization options structure
% Edited by Yi-Hsuan Chen 08/30/2023

%%========================= Optimization settings =========================
options.max_iterations  = 10;
options.converge_tol    = 1e-3;
options.clean_tol       = 1e-6;
options.degL1           = 2;
options.degL0           = 4;
options.degLe           = 2;
options.degV            = 4;
options.xdim            = length(x0(1));
options.rho0            = 0.1*ones(length(ts),1);
options.gX0             = [4 0 0; 0 4 0; 0 0 8];

%%============================ Initialize rho =============================
c           = 3;
rhot        = exp(c*(ts-max(ts))/(max(ts)-min(ts)));
rho0pp      = interp1(ts,rhot,'linear','pp');

%%================ convert piecewise polynomial to msspoly ================
N           = length(ts);
t           = msspoly('t',1);
x           = msspoly('x',options.xdim);
u0i         = pp2msspoly(t,spline(ts,u0(ts)));
x0i         = pp2msspoly(t,spline(ts,x0(ts)));
Spp         = interp1(ts,reshape(permute(S0(ts),[3 1 2]),N,options.xdim*options.xdim),'linear','pp');
S0i         = pp2msspoly(t,Spp);
% rhoi        = pp2msspoly(t,rhopp);

%===================== compute (V0,V0dot,rho,rhodot) ======================
intervals   = cell(1,N-1);
V0          = cell(1,N-1);
V0dot       = cell(1,N-1);
S0s         = cell(1,N-1);
% rho0         = cell(1,N-1);
% rhodot      = cell(1,N-1);
% Phi         = cell(1,N-1);
ui          = cell(1,N-1);
% rho         = ppval(rhopp,ts);

for i = 1:N-1  
    intervals{i}    = [ts(i) ts(i+1)]-ts(i);
%     Phi{i}          = zeros(length(x),length(x));   % initialize Phi
    S0s{i}          = reshape(S0i(:,i),options.xdim,options.xdim);
    xbari           = x-x0i(:,i);
    ui{i}           = u0i(:,i)+K(ts(i))*xbari;
    % Use matrix 'S' for initial Lyapunov guess
%     V0{i}           = xbari'*S0s{i}*xbari + x'*Phi{i}*x;
    V0{i}           = xbari'*S0s{i}*xbari;
    V0dot{i}        = diff(V0{i},t);
%     rho0{i}         = rhoi(:,i);
end
Vtraj0.V0       = V0;
% Vtraj0.V0dot    = V0dot;
Vtraj0.S0s      = S0s;
intvs           = intervals;
end