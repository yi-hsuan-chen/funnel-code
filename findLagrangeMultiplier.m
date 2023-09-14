%% Back up script for findL function
c       = 3;
rhot    = exp(c*(tGrid-max(tGrid))/(max(tGrid)-min(tGrid)));
rhopp   = interp1(tGrid,rhot,'linear','pp');
% convert piecewise polynomial to msspoly 
N       = length(tGrid);
t       = msspoly('t',1);
x       = msspoly('x',xdim);
u0i     = pp2msspoly(t,spline(tGrid,u0(tGrid)));
x0i     = pp2msspoly(t,spline(tGrid,x0(tGrid)));
SGrid   = S(tGrid);
Spp     = interp1(tGrid,reshape(permute(SGrid,[3 1 2]),N,xdim*xdim),'linear','pp');
S0i     = pp2msspoly(t,Spp);
rhoi    = pp2msspoly(t,rhopp);

% compute xdot{i}, V{i}, and rho{i}
intervals   = cell(1,N-1);
xdot        = cell(1,N-1);
V           = cell(1,N-1);
rho         = cell(1,N-1);
Phi         = cell(1,N-1);

for i = 1:N-1
t_now           = tGrid(i);
t_next          = tGrid(i+1);       
Phi{i}          = zeros(length(x),length(x));
intervals{i}    = [t_now t_next]-t_now;
S0s             = reshape(S0i(:,i),xdim,xdim);
xbari           = x-x0i(:,i);
xdot{i}         = f_polysym(t,x,u0i(:,i)+K(t_now)*xbari);
% Use matrix 'S' for initial Lyapunov guess
V{i}            = xbari'*S0s*xbari;
rho{i}          = rhoi(:,i);
end

k       = 2;
Vs      = repmat(V{1},N-1,k);
rhos    = repmat(V{1},N-1,k);
Vdots   = repmat(V{1},N-1,k);
rhodots = repmat(V{1},N-1,k);

for i = 1:N-1
    tt          = linspace(intervals{i}(1),intervals{i}(2),k);
    tts(i,:)    = tt;
    Vs(i,:)     = msubs(V{i},t,tt);
    
    % Remember(need to think about how to add) Phik!!]
%     Phidoti = (Phi{i+1}-Phi{i})/(tGrid(i+1)-tGrid(i));
    Vdotsi      = diff(V{i},t)+diff(V{i},x)*xdot{i};
    Vdots(i,:)  = msubs(Vdotsi,t,tt);
    rhos(i,:)   = msubs(rho{i},t,tt);
    rhodots(i,:)= msubs(diff(rho{i},t),t,tt);
end

%%
M   = length(Vs(:));
disp('Step 1: Searching for multipliers...')
L1f = cell(1,M);

for i = 1:M

    % Initialize SOS program
    prog    = spotsosprog;
    prog    = prog.withIndeterminate(x);
    Vi      = clean(Vs(i),options.clean_tol);
    Vdoti   = clean(Vdots(i),options.clean_tol);  
    rhoi    = rhos(i);
    rhodoti = rhodots(i);
    
    % create multiplier and slack variable
    [prog,L1]   = prog.newFreePoly(monomials(x,0:options.degL1));
    [prog,slack]= prog.newPos(1);
    
    conDeriv    = -slack*(x'*x)^4 - Vdoti + rhodoti + L1*(Vi-rhoi);
    prog        = prog.withSOS(conDeriv);
    
    % solve problem
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    costfun                 = -slack;
    sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
    
    if ~strcmp(sol.status,'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
        disp('Primal infeasible')
    end
    L1f{i}      = sol.eval(L1);
    slacks(i)   = double(sol.eval(slack));
end