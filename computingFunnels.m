clc; clear; close all;

%% Step 1: Trajectory Generation
traj_data   = load('trajectories.mat');
traj_lib    = traj_data.solutions;
traj_lib    = squeeze(num2cell(traj_lib,[1 2]));
traj_grid   = traj_data.solutionsGrid;
traj_grid   = squeeze(num2cell(traj_grid,[1 2]));
nSample     = size(traj_grid{1},2);
xdim        = 3;
udim        = 1;

fig(1) = figure(); ax(1) = gca; hold on; grid on;
for i = 1:length(traj_lib)
    traj_temp       = traj_lib{i}';
    traj_grid_temp  = traj_grid{i}';
    t_nom           = traj_temp(:,1);
    x_nom           = traj_temp(:,2);
    y_nom           = traj_temp(:,3);
    th_nom          = traj_temp(:,4);
    u_nom           = traj_temp(:,end);
    tGrid           = traj_grid_temp(:,1);
    xGrid           = traj_grid_temp(:,2);
    yGrid           = traj_grid_temp(:,3);
    thGrid          = traj_grid_temp(:,4);
    uGrid           = traj_grid_temp(:,end);
    % Plot the entire trajectory
    plot(x_nom,y_nom,'r-','LineWidth',3);
    % Plot the grid points:
    plot(xGrid, yGrid, 'ko','MarkerSize',3,'LineWidth',3);
    % Plot the start and end points:
    plot(x_nom([1,end]), y_nom([1,end]),'ks','MarkerSize',12,'LineWidth',3);
end

%% Step 2: Dynamics Setup
% TODO: 
% 0) Fit (xo,uo) to piecewise polynomials
X_nom       = [x_nom y_nom th_nom];
Xpp         = spline(t_nom,X_nom');
Upp         = spline(t_nom,u_nom);

% 1) Linearization the system around the nominal trajectory (xo,uo)
paras.v             = 1;
[state,f,g,fh,gh]   = dubins_dynamics(paras);
x0                  = @(t) ppval(Xpp,t);
u0                  = @(t) ppval(Upp,t);
[A,B]               = linearizeSys(state,f,g,x0,u0);

% 2) Compute TVLQR by solving matrix Riccati equation (backward in time!)
Q           = @(t) diag([10 10 1]);
R           = @(t) 1;
S0          = diag([2 2 4]);
tspan       = [0 t_nom(end)];
[ts,Ss]     = tvlqr_riccati(tspan,A,B,Q,R,S0);
% Spp         = spline(ts,Ss);
S           = @(t) ppval(spline(ts,Ss),t);
K           = @(t) -R(t)\B(t)'*S(t);

%% 3) Taylor expanded the CL sytem around (xo,uo) to get polynomial dynamics
fcl         = @(t,x) fh(x) + gh(x)*(u0(t)+K(t)*(x-x0(t)));
checkTVLQRctrlperf(fcl,tspan,X_nom);
order       = 3;
% fcl_poly    = taylorExpansion(state,f,g,K,x0,order);
% [~,x_app]   = ode45(fcl_poly,tspan,[0;0;0]);
f_polysym   = taylorApproxSyms(fh,gh,order);
[~,x_app]   = ode45(@(t,x) f_polysym(t,x,u0(t)+K(t)*(x-x0(t))),tspan,[0;0;0]);
plot(x_app(:,1),x_app(:,2),'k:','LineWidth',3);

%% Step 3: Implement Algorithm 1

%
% 3.1) Initialization: V0,V0dot,rho,rhodot,Phi,ui
% ui is the nominal control expressed in msspoly
[Vtraj0,rho0pp,Phi,ui,intvs,t,x,options] = InitializeOptim(tGrid,x0,u0,K,S);

% [L1f,slackf]    = findL(tGrid,f_polysym,x0,u0,K,S,options);

% % 3.2) Perform bilinear alternation search
V               = Vtraj0.V0;
rhopp           = rho0pp;
k               = 1;

[L1,Le,slacks]  = findMultipliers(V,rhopp,tGrid,f_polysym,Phi,options,Vtraj0,ui,intvs,t,x,k);

% for iter = 1:1
    % First step: Fix V and rho, search for L
%     [L1,slacks]  = findL1(V,rhopp,tGrid,f_polysym,Phi,options,Vtraj0,ui,intvs,t,x,k);
    % Second step: Fix L, search for V and rho (and Phi)
    % rho should be returned as a piecewise polynomial NOT double array
    %     Phiold  = Phi; 
    %     [V,rhopp,Phi]   = findVRho(V,L1,tGrid,f_polysym,Phi,options,Vtraj0,ui,intvs,t,x)
% end

%% Step 2: Fix L, search for rho and V (Phi)
N       = length(tGrid)-1;
% Initialize program
prog    = spotsosprog;
prog    = prog.withIndeterminate(x);

% Get (V0,V0dot,S0) at each time step -- remove time-dependency
V0s         = repmat(Vtraj0.V0{1},N,k);
V0dots      = repmat(Vtraj0.V0{1},N,k);
xdots       = cell(N,k);
S0s         = cell(N,k);
for i = 1:N
    tt          = linspace(intvs{i}(1),intvs{i}(2),k);
    V0s(i,:)    = msubs(Vtraj0.V0{i},t,tt);
    V0dots(i,:) = msubs(Vtraj0.V0dot{i},t,tt);
    xdot_temp   = msubs(f_polysym(t,x,ui{i}),t,tt);
    xdots(i,:)  = {xdot_temp(:,1:k)};
    S0s_temp    = msubs(reshape(Vtraj0.S0s{i},options.xdim*options.xdim,1),t,tt);
    S0s(i,:)    = {full(reshape(S0s_temp,options.xdim,options.xdim))};
end

%%
M   = length(V0s(:));
disp('Step 2: Searching for V(Phid) and rho...')
% Initialize rho and PhiM
[prog,rho]          = prog.newPos(M);

% check this constr...
prog                = prog.withPos(rho(1)-ppval(rhopp,tGrid(1)));
rhodot              = msspoly(zeros(M,1));
PhidM               = cell(M,1);

for i = 1:M
    % Make V positive definite
    [prog,PhidM{i}] = prog.newSym(length(x));
    prog            = prog.withPSD(S0s{i}+PhidM{i});
end
[prog,tau]          = prog.newPos(1);
VV{1}               = V0s(1) + x'*PhidM{1}*x;

% check this constr...
prog                = prog.withSOS(rho(1)-VV{1}-tau*(x'*S0*x));

%
dt      = intvs{1}(2)-intvs{1}(1);
for i = 1:M
    Vdoti           = diff(V0s(i),x)*xdots{i} + V0dots(i) + 2*x'*PhidM{i}*xdots{i};
    Vi              = V0s(i) + x'*PhidM{i}*x;
    
    % Use forward/backward for boundary and centered diff for the rest
    if i <= 2
        rhodot(i)       = (-3*rho(i)+4*rho(i+1)-rho(i+2))/(2*dt);
    elseif i >= M-1
        rhodot(i)       = (3*rho(i)-4*rho(i-1)+rho(i-2))/(2*dt);
    else
        rhodot(i)       = (-rho(i+2)+8*rho(i+1)-8*rho(i-1)+rho(i-2))/(12*dt);
    end
    
    % clean stuff
    L1{i}           = clean(L1{i},options.clean_tol);
    Vi              = clean(Vi,options.clean_tol);
    Vdoti           = clean(Vdoti,options.clean_tol);
    
    conDeriv        = -Vdoti + rhodot(i) - L1{i}*(Vi-rho(i));
    prog            = prog.withSOS(conDeriv);  
    
    % new constr...
    [prog,Le{i}]    = prog.newFreePoly(monomials(x,0:1));
    prog            = prog.withSOS(Le{i});
    conderiv        = 1-x'*PhidM{i}*x-Le{i}*(rho(i)-Vi);
    prog            = prog.withSOS(conderiv);
    
    
    % Solve SOS program
    options_spot            = spot_sdp_default_options();
    options_spot.verbose    = 0;
    costfun                 = rho(i)/1000;
    sol                     = prog.minimize(costfun,@spot_mosek,options_spot);

    if ~strcmp(sol.status,'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
        disp('Primal infeasible')
    end
    rho_opt(i)  = double(sol.eval(rho(i)))
end

% Solve SOS program
% options_spot            = spot_sdp_default_options();
% options_spot.verbose    = 0;
% costfun                 = sum(rho)/10000000;
% sol                     = prog.minimize(costfun,@spot_mosek,options_spot);
% 
% if ~strcmp(sol.status,'STATUS_PRIMAL_AND_DUAL_FEASIBLE')
%     disp('Primal infeasible')
% end
% rho_opt  = double(sol.eval(rho))
% for i = 1:M-1
%     Phi_opt{i}  = double(sol.eval(PhidM{i}));
% end



function checkTVLQRctrlperf(fcl,tspan,X_nom)
% Note: Check control performance of TVLQR for nonlinear sys
x0          = [0;0;0];
[~,x_n]     = ode45(fcl,tspan,x0);
x_nom       = X_nom(:,1);
y_nom       = X_nom(:,2);
plot(x_nom,y_nom,'r-',x_n(:,1),x_n(:,2),'b--','LineWidth',3); grid on; 
end