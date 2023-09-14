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

[L1,Le,Sks]     = findLs(V,rhopp,tGrid,f_polysym,options,ui,intvs,t,x,k);


function checkTVLQRctrlperf(fcl,tspan,X_nom)
% Note: Check control performance of TVLQR for nonlinear sys
x0          = [0;0;0];
[~,x_n]     = ode45(fcl,tspan,x0);
x_nom       = X_nom(:,1);
y_nom       = X_nom(:,2);
plot(x_nom,y_nom,'r-',x_n(:,1),x_n(:,2),'b--','LineWidth',3); grid on; 
end