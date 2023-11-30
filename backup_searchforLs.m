function [L1s,LEs,Sks,slacks] = searchforLs(V,rho,ts,xdots,options,x)
% Goal: Minimize vol(E(tk)) by searching for (L1,LE,L0,Sk)
% 
% Input:
%   V       -- Lyapunov fcn in 1-by-(N-1) cell with msspoly (independent of t)
%   rho     -- upper bound of Lyapunov fcn (independent of t)
%   ts      -- time samples
%   xdots   -- Taylor expanded CL system in error coordinate (xbardot in the text)
%   options -- optimization options structure
% % %   Vtraj0  -- structure including V0, V0k(dV0/dt), S0s expressed in 1-by-(N-1) cells with msspoly
% % %   rho0    -- struture including rho and rhodot expressed in 1-by-(N-1) cells with msspoly
% % %   ui      -- tvlqr CL control input expressed in 1-by-(N-1) cells with msspoly
% % %   intvs   -- time intervals expressed in 1-by-(N-1) cell
% % %   t       -- time expressed in 1-by-1 msspoly 
%   x       -- state expressed in n-by-1 msspoly
% % %   k       -- number of time segments
% Output:
%   L1s     -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
% %   LEs     -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
% %   SKs     -- PSD matrices associated with the volume of the ellipsoid
%   (slacks)-- slack (double) variables for debugging
% Edited by Yi-Hsuan Chen 10/02/2023

%%=========================== Initialization ===========================
N       = length(ts);
n       = N-1;
Vdot    = cell(1,n);
rho_    = zeros(1,n);
dt      = ts(2)-ts(1);         % step size for numerical differentiation

for i = 1:length(rho)
    rho_(i) = double(rho{i});
end
rhodot  = gradient(rho_);

for i = 1:n
    xdoti       = xdots{i};
    dVdx        = diff(V{i},x);
    switch i
        case 1
            % use FORWARD difference here for the first point
            dVdt        = (-3*V{i}+4*V{i+1}-V{i+2})/(2*dt);
        case n
            % use BACKWARD difference here for the last point
            dVdt        = (3*V{i}-4*V{i-1}+V{i-2})/(2*dt);
        otherwise
            % use CENTRAL difference
            dVdt        = (V{i+1}-V{i-1})/(2*dt);
    end
    Vdot{i}     = dVdx*xdoti + dVdt;    
end

L1s     = cell(1,n);
LEs     = cell(1,n);
% SKs     = cell(1,n);
% slacks  = zeros(1,n);
fprintf('Step 1: Searching for multipliers...')

% % Initialize SOS program
% prog        = spotsosprog;
% prog        = prog.withIndeterminate(x);

% [prog,L0]   = prog.newSOSPoly(monomials(x,0:options.degL0));
%         prog        = prog.withSOS(L0);
% constr3     = rho_(1)-V{1}-L0*(1-x'*options.gX0*x);
% prog        = prog.withSOS(constr3);

% create multipliers and slack variable
% [prog,slack]= prog.newFree(n);
% [prog,L1]   = prog.newFreePoly(monomials(x,0:options.degL1),n);

for i = 1:n
%     i
    % Initialize SOS program
    prog        = spotsosprog;
    prog        = prog.withIndeterminate(x);

    % clean stuff
    Vi          = clean(V{i},options.clean_tol);
    Vdoti       = clean(Vdot{i},options.clean_tol);
    rhoi        = rho_(i);
    rhodoti     = rhodot(i);
    
    % create multipliers and slack variable
%     [prog,gamma]= prog.newFree(1);
    [prog,L1]   = prog.newFreePoly(monomials(x,0:options.degL1));

% %     constr1     = slack(i) -Vdoti + rhodoti - L1(i)*(Vi-rhoi);
%     constr1     = gamma -Vdoti + rhodoti - L1*(Vi-rhoi);
    constr1     = -Vdoti + rhodoti + L1*(Vi-rhoi);
    prog        = prog.withSOS(constr1);
    
%     Create a PSD matrix Sk
    [prog,Sk]   = prog.newPSD(options.xdim);
%     prog        = prog.withPSD(Sk);
    [prog,Le]   = prog.newSOSPoly(monomials(x,0:options.degLe));
%     prog        = prog.withSOS(Le);
    constr2     = 1-x'*Sk*x-Le*(rhoi-Vi);
    prog        = prog.withSOS(constr2);

%     if i==1
%         [prog,L0]   = prog.newSOSPoly(monomials(x,0:options.degL0));
% %         prog        = prog.withSOS(L0);
% %         [prog,L0]   = prog.newPos(1);
%         constr3     = rho_(1)-V{1};
%         prog        = prog.withSOS(constr3);
%     end
% end

% solve problem
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    
    %%... check the correctness of maxdet function
    [prog,obj]              = maxdet(prog,Sk);
    costfun                 = -obj;
%     costfun                 = gamma;
%     costfun                 = msspoly(0);
    try
        sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
    catch
        % failed
        fprintf('failed \n');
        return;
    end
    solved      = sol.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');
    % Parse
    if solved
        fprintf('solved \n');
        clean_eps   = 1e-8;
        L1s{i}      = clean(sol.eval(L1),clean_eps);
        LEs{i}      = clean(sol.eval(Le),clean_eps);
        Sks{i}      = clean(double(sol.eval(Sk)),clean_eps);
%         slacks(i)   = clean(double(sol.eval(gamma)),clean_eps);
        clear prog;
    else
        % failed
        fprintf('failed \n');
        return;
    end    
end
    slacks      = 0;
end