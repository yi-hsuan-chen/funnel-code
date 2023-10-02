%%=========================== Initialization ===========================
N       = length(TGrid);
rhoi    = pp2msspoly(t,rhopp);
Vs      = repmat(V{1},N-1,k);
Vdots   = repmat(V{1},N-1,k);
rhos    = repmat(V{1},N-1,k);
rhodots = repmat(V{1},N-1,k);


for i = 1:N-1
    tt          = linspace(intvs{i}(1),intvs{i}(2),k);
    Vs(i,:)     = subss(V{i},t,tt);
    xdoti       = polydyn(t,x,ui{i});
    Vdotsi      = diff(V{i},x)*xdoti + diff(V{i},t);
    Vdots(i,:)  = subss(Vdotsi,t,tt);
    rhos(i,:)   = subss(rhoi(:,i),t,tt);
    rhodots(i,:)= subss(diff(rhoi(:,i),t),t,tt);
end

L1s     = cell(1,N-1);
Les     = cell(1,N-1);
Sks     = cell(1,N-1);
slacks  = zeros(1,N-1);
for i = 1:N-1
    disp('Step 1: Searching for multipliers...')
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
    constr1     = -Vdoti + rhodoti - L1*(Vi-rhoi);
    prog        = prog.withSOS(constr1);
    
    % Create a PSD matrix Sk
    [prog,Sk]   = prog.newPSD(options.xdim);
    [prog,Le]   = prog.newFreePoly(monomials(x,0:options.degLe));
    prog        = prog.withSOS(Le);
    constr2     = 1-x'*Sk*x-Le*(rhoi-Vi);
    prog        = prog.withSOS(constr2);
    
%     [prog,L0]   = prog.newFreePoly(monomials(x,0:options.degL0));
%     prog        = prog.withSOS(L0);
%     constr3     = rhos(1)-Vs(1)-L0*(1-x'*options.gX0*x);
%     prog        = prog.withSOS(constr3);

    % solve problem
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    
    %%... check the correctness of maxdet function
    [prog,obj]              = maxdet(prog,Sk);
    costfun                 = -obj;
    sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
    
    if ~sol.isPrimalFeasible || ~sol.isDualFeasible
        disp('Infeasibility detected')
    end
    L1s{i}      = sol.eval(L1);
    Les{i}      = sol.eval(Le);
    Sks{i}      = double(sol.eval(Sk));
end