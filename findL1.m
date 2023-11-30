function [L1s,slacks] = findL1(V,rho,ts,xdots,Vtraj0,Phi,options,x)
% Goal: Minimize vol(E(tk)) by searching for (L1,LE,L0,Sk)
% 
% Input:
%   V       -- Lyapunov fcn in 1-by-(N-1) cell with msspoly (independent of t)
%   rho     -- upper bound of Lyapunov fcn (independent of t)
%   ts      -- time samples
%   xdots   -- Taylor expanded CL system in error coordinate (xbardot in the text)
%   Vtraj0  == structure including V0, V0k(dV0/dt)
%   Phi     == 
%   options -- optimization options structure
%   x       -- state expressed in n-by-1 msspoly
% % %   Vtraj0  -- structure including V0, V0k(dV0/dt), S0s expressed in 1-by-(N-1) cells with msspoly
% % %   rho0    -- struture including rho and rhodot expressed in 1-by-(N-1) cells with msspoly
% % %   ui      -- tvlqr CL control input expressed in 1-by-(N-1) cells with msspoly
% % %   intvs   -- time intervals expressed in 1-by-(N-1) cell
% % %   t       -- time expressed in 1-by-1 msspoly 
% % %   k       -- number of time segments
%
% Output:
%   L1s     -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   LEs     -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   SKs     -- PSD matrices associated with the volume of the ellipsoid
%   (slacks)-- slack (double) variables for debugging
% Edited by Yi-Hsuan Chen 10/02/2023

%%=========================== Initialization ===========================
N       = length(ts);
n       = N-1;
Vdot    = cell(1,n);
rho_    = zeros(1,n);
dt      = ts(2)-ts(1);          % step size for numerical differentiation

for i = 1:length(rho)
    rho_(i) = double(rho{i});   % convert rho (msspoly) into rho_ (double)
end
rhodot  = gradient(rho_);

for i = 1:n
    xdoti       = xdots{i};
    V0doti      = Vtraj0.V0dot{i};
    Phidoti     = diff_wrt_t(i,Phi,dt,n);
    Vdot{i}     = V0doti + 2*x'*Phi{i}*xdoti + x'*Phidoti*x;
end
L1s     = cell(1,n);
slacks  = zeros(1,n);

fprintf('Step 1: Searching for multipliers...\n')
for i = 1:n
    fprintf('iter: %d...',i);
    
    % Initialize SOS program
    prog        = spotsosprog;
    prog        = prog.withIndeterminate(x);

    % clean stuff
    Vi          = clean(V{i},options.clean_tol);
    Vdoti       = clean(Vdot{i},options.clean_tol);
    rhoi        = rho_(i);
    rhodoti     = rhodot(i);
    
    % create multipliers and slack variable
    [prog,gamma]= prog.newFree(1);
    [prog,L1]   = prog.newFreePoly(monomials(x,0:options.degL1));
    constr1     = gamma -Vdoti + rhodoti - L1*(Vi-rhoi);
    prog        = prog.withSOS(constr1);

    % solve problem
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    costfun                 = gamma;
    try
        sol                 = prog.minimize(costfun,@spot_mosek,options_spot);
    catch
        % failed
        fprintf('failed \n');
        return;
    end
    solved      = sol.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');
    
    % Parse
    if solved
        fprintf('solved \n');
        L1s{i}      = clean(sol.eval(L1),options.clean_tol);
        slacks(i)   = clean(double(sol.eval(gamma)),options.clean_tol);
        clear prog;
    else
        % failed
        fprintf('failed \n');
        return;
    end    
end
end