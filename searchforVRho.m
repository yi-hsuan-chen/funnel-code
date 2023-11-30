function [VV,Rho,Phis] = searchforVRho(L1,Le,rho_old,gamma,ts,xdots,Vtraj0,Phi_old,options,x)
% Goal: Minimize vol(E(tk)) by searching for (V,rho,L0,Phis)
% 
% Input:
%   L1      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   Le      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   rho_old -- upper bound of Lyapunov fcn (scalar?)
%   ts      -- time samples
%   xdots   -- Taylor expanded CL system in error coordinate (xbardot in the text)
%   options -- optimization options structure
%   x       -- state expressed in n-by-1 msspoly
%
% Output:
%   VV      -- Lyapunov fcn in 1-by-(N-1) cell with msspoly (independent of t) 
%   Rho     -- upper bound of Lyapunov fcn (scalar?)
%   SKs     -- PSD matrices associated with the volume of the ellipsoid
% Edited by Yi-Hsuan Chen 10/04/2023

dts         = diff(ts);
N           = length(ts);
n           = N-1;
Vs          = cell(1,n);
rhos        = cell(1,n);
Vdots       = cell(1,n);
rhodots     = cell(1,n);
Sk          = cell(1,n);

fprintf('Step 2: Searching for V and rho...');
% Initialize SOS program
prog        = spotsosprog;
prog        = prog.withIndeterminate(x);
for k = 1:n
    [prog,rhos{k}]  = prog.newPos(1);
%     prog            = prog.withPos(rhos{k}-rho_old{k});
%     [prog,Vs{k}]    = prog.newFreePoly(monomials(x,0:options.degV));
%     [prog,Sk{k}]    = prog.newPSD(options.xdim);
    [prog,Phi{k}]   = prog.newSym(options.xdim);
    Sk              = Vtraj0.S0s{k}+Phik;
    prog            = prog.withPSD(Sk);
%     prog            = prog.withPSD(Sk{k}-options.delta_s*eye(options.xdim));
    V{k}        = Vtraj0.V0{k} + x'*Phi{k}*x;
    V0dotk      = Vtraj0.V0dot{k};
    Phidotk     = diff_wrt_t(k,Phi,dt,n);
    Vdot{k}     = V0doti + 2*x'*Phi{k}*xdotk + x'*Phidotk*x;
end
prog        = prog.withEqs(rhos{1}-rho_old{1});
% % prog        = prog.withEqs(rhos{n}-1);

% [prog,E_eps]= prog.newFree(1);
[prog,L0]   = prog.newSOSPoly(monomials(x,0:options.degL0));
constr      = rhos{1}-Vs{1}-L0*(1-x'*options.gX0*x);
prog        = prog.withSOS(constr);

costfun     = 0; 
for k = 1:n
    dVdx        = diff(Vs{k},x);
    dVdt        = diff_wrt_t(k,Vs,dts(k),n);
    Vdots{k}    = dVdx*xdots{k} + dVdt; 
    rhodots{k}  = diff_wrt_t(k,rhos,dts(k),n);

    % constraints
    constr1     = -Vdots{k}+rhodots{k}-L1{k}*(Vs{k}-rhos{k});
    prog        = prog.withSOS(constr1);

    constr2     = 1-x'*Sk{k}*x-Le{k}*(rhos{k}-Vs{k});
    prog        = prog.withSOS(constr2);
    
    [prog,obj]  = maxdet(prog,Sk{k});
    costfun     = costfun-obj;
end

% solve problem
options_spot            = spot_sdp_default_options();
options_spot.verbose    = 0;
try
    sol                 = prog.minimize(costfun,@spot_mosek,options_spot);
catch
    % failed
    fprintf('failed \n');
    return;
end
solved                  = sol.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');
if solved
    for i = 1:n
        VV{i}       = clean(sol.eval(Vs{i}),options.clean_tol);
        Rho{i}      = clean(double(sol.eval(rhos{i})),options.clean_tol);
        Sks{i}      = clean(double(sol.eval(Sk{i})),options.clean_tol);
    end
    clear prog;
else
    % failed
    fprintf('failed \n');
    return;
end
end