function [VV,Rho,Phi_] = findVRho(L1,rho_old,gamma,ts,xdots,Vtraj0,Phi_old,options,x)
% Goal: Minimize vol(E(tk)) by searching for (V,rho,L0,Phis)
% 
% Input:
%   L1      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
% % %   Le      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   rho_old -- upper bound of Lyapunov fcn (scalar?)
%   ts      -- time samples
%   xdots   -- Taylor expanded CL system in error coordinate (xbardot in the text)
%   options -- optimization options structure
%   x       -- state expressed in n-by-1 msspoly
%
% Output:
%   VV      -- Lyapunov fcn in 1-by-(N-1) cell with msspoly (independent of t) 
%   Rho     -- upper bound of Lyapunov fcn (scalar?)
%   Phi_    -- sym. matrices associated with the volume of the ellipsoid
% Edited by Yi-Hsuan Chen 10/04/2023

dts         = diff(ts);
N           = length(ts);
n           = N-1;
V           = cell(1,n);
rho         = cell(1,n);
Vdot        = cell(1,n);
rhodot      = cell(1,n);
PhidM       = cell(1,n);
Sk          = cell(1,n);

fprintf('Step 2: Searching for V and rho...');
% % Initialize SOS program
% prog            = spotsosprog;
% prog            = prog.withIndeterminate(x);
% [prog,PhidM{1}] = prog.newSym(length(x));
% V{1}            = Vtraj0.V0{1} + x'*PhidM{1}*x;

% for i = 1:n
%     [prog,rho{i}] = prog.newPos(1);
% end
% prog            = prog.withPos(rho{1}-rho_old{1});

% [prog,L0]       = prog.newSOSPoly(monomials(x,0:options.degL0));
[prog,L0]       = prog.newPos(1);
constr          = rho{1}-V{1}-L0*(1-x'*options.gX0*x);
% prog            = prog.withSOS(constr);

for k = 2:n
    [prog,PhidM{k}] = prog.newPSD(options.xdim);
    Sk{k}           = Vtraj0.S0s{k}+PhidM{k};
%     prog            = prog.withPSD(Sk{k});
end

costfun     = 0; 
% C           = [eye(2),zeros(2,1)];
for k = 1:n
    V{k}        = Vtraj0.V0{k} + x'*PhidM{k}*x;
    Phidotk     = diff_wrt_t(k,PhidM,dts(k),n);
    Vdot{k}     = Vtraj0.V0dot{k} + 2*x'*PhidM{k}*xdots{k} + x'*Phidotk*x;
    rhodot{k}   = diff_wrt_t(k,rho,dts(k),n);

    % clean stuff
    L1{k}       = clean(L1{k},options.clean_tol);
    V{k}        = clean(V{k},options.clean_tol);
    Vdot{k}     = clean(Vdot{k},options.clean_tol);

    % constraints
    constr1     = -Vdot{k}+rhodot{k}-L1{k}*(V{k}-rho{k});
    prog        = prog.withSOS(constr1);
    
%     costfun     = costfun + rho{k};

    [prog,obj]  = maxdet(prog,Sk{k});
    costfun     = costfun-obj;
%     P0k         = (Vtraj0.S0s{k}+Phi_old{k})/double(rho{k});
%     Pk          = (Vtraj0.S0s{k}+PhidM{k})/double(rho{k});
%     costfun     = costfun-trace(C'*inv(C*inv(P0k)*C')*C*P0k*Pk*inv(P0k));
end
% solve problem
options_spot                    = spot_sdp_default_options();
options_spot.verbose            = 0;
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
        VV{i}       = clean(sol.eval(V{i}),options.clean_tol);
        Rho{i}      = clean(double(sol.eval(rho{i})),options.clean_tol);
        Phi_{i}      = clean(double(sol.eval(PhidM{i})),options.clean_tol);
    end
    clear prog;
else
    % failed
    fprintf('failed \n');
    return;
end
end