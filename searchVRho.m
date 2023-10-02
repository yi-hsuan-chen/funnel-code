% Search for V and rho
dts         = diff(TGrid);
N           = length(TGrid);
disp('Step 2: Searching for V and rho...');
% Initialize SOS program
prog        = spotsosprog;
prog        = prog.withIndeterminate(x);
[prog,rho]  = prog.newPos(N);
prog        = prog.withPos(rho(1)-ppval(rho0pp,0));

for k = 1:N
    [prog,Vs{k}]    = prog.newFreePoly(monomials(x,0:options.degV));
end
costfun
for k = 1:N-1
    Vdots{k}    = (Vs{k+1}-Vs{k})/dts(k);
    rhodots{k}  = (rho(k+1)-rho(k))/dts(k);

    % create multipliers
    constr1     = -Vdots{k} + rhodots{k} - L1{k}*(Vs{k+1}-rho(k+1));
    prog        = prog.withSOS(constr1);

    % Create a PSD matrix Sk
    [prog,Sk]   = prog.newPSD(options.xdim);
    constr2     = 1-x'*Sk*x-Le{k}*(rho(k+1)-Vs{k+1});
    prog        = prog.withSOS(constr2);
    %
    %             [prog,L0]   = prog.newFreePoly(monomials(x,0:options.degL0));
    %             prog        = prog.withSOS(L0);
    %             constr3     = rho(1)-Vs{1}-L0*(1-x'*options.gX0*x);
    %             prog        = prog.withSOS(constr3);

    [prog,obj]  = maxdet(prog,Sk);
    costfun     = costfun-obj;
end
% solve problem
options_spot            = spot_sdp_default_options();
options_spot.verbose    = 0;
sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
sol.status
if ~sol.isPrimalFeasible || ~sol.isDualFeasible
    disp('Infeasibility detected')
end
VV      = sol.eval(Vs);
Rho     = sol.eval(rho);
Skks    = double(sol.eval(Sk));