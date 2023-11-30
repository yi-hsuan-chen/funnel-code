% Search for V and rho
dts         = diff(TGrid);
N           = length(TGrid);
n           = N-1;
Vs          = cell(1,n);
rhos        = cell(1,n);
Vdots       = cell(1,n);
rhodots     = cell(1,n);
Skk         = cell(1,n);

fprintf('Step 2: Searching for V and rho...');
% Initialize SOS program
prog        = spotsosprog;
prog        = prog.withIndeterminate(x);
for k = 1:n
    [prog,rhos{k}]  = prog.newPos(1);
    prog            = prog.withPos(rhos{k}-rho{k});

    [prog,Vs{k}]    = prog.newFreePoly(monomials(x,0:options.degV));
    [prog,Skk{k}]   = prog.newPSD(options.xdim);
%     prog            = prog.withPSD(Skk{k}-options.delta_s*eye(options.xdim));
end
% prog        = prog.withEqs(rhos{1}-rho{1});
% prog        = prog.withEqs(rhos{n}-1);

[prog,L0]   = prog.newSOSPoly(monomials(x,0:options.degL0));
constr3     = rhos{1}-Vs{1}-L0*(1-x'*options.gX0*x);
prog        = prog.withSOS(constr3);

costfun     = 0; 
for k = 1:n
    dVdx        = diff(Vs{k},x);
    dVdt        = diff_wrt_t(k,Vs,dts(k),n);
    Vdots{k}    = dVdx*xdots{k} + dVdt; 
    rhodots{k}  = diff_wrt_t(k,rhos,dts(k),n);

    % constraints
    constr1     = -Vdots{k}+rhodots{k}-L1{k}*(Vs{k}-rhos{k});
    prog        = prog.withSOS(constr1);

    constr2     = 1-x'*Skk{k}*x-Le{k}*(rhos{k}-Vs{k});
    prog        = prog.withSOS(constr2);
    
    [prog,obj]  = maxdet(prog,Skk{k});
    costfun     = costfun-obj;
end

% solve problem
options_spot            = spot_sdp_default_options();
options_spot.verbose    = 0;
sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
solved                  = sol.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');
if solved
    for i = 1:n
        VVs{i}      = sol.eval(Vs{i});
        Rhos{i}     = double(sol.eval(rhos{i}));
        Skks{i}     = double(sol.eval(Skk{i}));
    end
    clear prog;
else
    % failed
    fprintf('failed \n');
    return;
end
