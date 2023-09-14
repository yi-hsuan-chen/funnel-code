function [L1s,Les,slacks] = findMultipliers(V,rhopp,ts,fpoly,Phi,options,Vtraj0,ui,intvs,t,x,k)
% 
% Input:
%   V       -- updated Lyapunov fcn (V0+x'*Phi*x) in 1-by-(N-1) cell with msspoly
%   rhopp   -- updated piecewise polynomial rho 
%   ts      -- time samples
%   fpoly   -- fpoly(t,x,u) Taylor expanded original system (not CL) 
%   Phi     -- Phi expressed in 1-by-(N-1) cell with n-by-n matrices 
%   options -- optimization options structure
%   Vtraj0  -- structure including V0, V0k(dV0/dt), S0s expressed in 1-by-(N-1) cells with msspoly
% % %   rho0    -- struture including rho and rhodot expressed in 1-by-(N-1) cells with msspoly
%   ui      -- tvlqr CL control input expressed in 1-by-(N-1) cells with msspoly
%   intvs   -- time intervals expressed in 1-by-(N-1) cell
%   t       -- time expressed in 1-by-1 msspoly 
%   x       -- state expressed in n-by-1 msspoly
%   k       -- number of time segments
% Output:
%   L1      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   Le      -- Lagrange multipliers in (N-1)-by-1 cell with msspoly
%   (slacks)-- slack (double) variables for debugging
% Edited by Yi-Hsuan Chen 08/30/2023

%%=========================== Initialization ===========================
N       = length(ts);
% V0      = Vtraj0.V0;
% V0dot   = Vtraj0.V0dot;

rhoi    = pp2msspoly(t,rhopp);
Vs      = repmat(V{1},N-1,k);
Vdots   = repmat(V{1},N-1,k);
rhos    = repmat(V{1},N-1,k);
rhodots = repmat(V{1},N-1,k);

for i = 1:N-1
    tt          = linspace(intvs{i}(1),intvs{i}(2),k);
    tts(i,:)    = tt;
    Vs(i,:)     = msubs(V{i},t,tt);
    xdoti       = fpoly(t,x,ui{i});
    Vdotsi      = diff(V{i},x)*xdoti + diff(V{i},t);
    Vdots(i,:)  = msubs(Vdotsi,t,tt);
    rhos(i,:)   = msubs(rhoi(:,i),t,tt);
    rhodots(i,:)= msubs(diff(rhoi(:,i),t),t,tt);
end

M       = length(Vs(:));
L1s     = cell(1,M);
Les     = cell(1,M);
disp('Step 1: Searching for multipliers...')

for i = 1:M
    % Initialize SOS program
    prog        = spotsosprog;
    prog        = prog.withIndeterminate(x);
    
    % clean stuff
    Vi          = clean(Vs(i),options.clean_tol);
    Vdoti       = clean(Vdots(i),options.clean_tol);  
    
    rhoi        = rhos(i);
    rhodoti     = rhodots(i);
    
    % create multipliers and slack variable
    [prog,L1]   = prog.newFreePoly(monomials(x,0:options.degL1));
    [prog,slack]= prog.newPos(1);
    constr1     = -slack*(x'*x)^4 - Vdoti + rhodoti + L1*(Vi-rhoi);
    prog        = prog.withSOS(constr1);
    
    [prog,Le]   = prog.newFreePoly(monomials(x,0:options.degLe));
    constr2     = 1-x'*Phi{i}*x-Le*(rhoi-Vi);
    prog        = prog.withSOS(constr2);
    
    % solve problem
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    costfun                 = -slack;
    sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
    
    if ~strcmp(sol.status,'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
        disp('Primal infeasible')
    end
    L1s{i}      = sol.eval(L1);
    Les{i}      = sol.eval(Le); 
    slacks(i)   = double(sol.eval(slack));
end

end