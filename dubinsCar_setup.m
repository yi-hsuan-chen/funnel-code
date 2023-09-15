clc; clear; close all;

%% Step 1: Trajectory Generation
load_trajectory_lib;

%% Step 2: Dynamics Setup
% 1) Fit (xo,uo) to piecewise polynomials
X_nom       = [x_nom y_nom th_nom];
Xpp         = spline(t_nom,X_nom');
Upp         = spline(t_nom,u_nom);
x0          = @(t) ppval(Xpp,t);
u0          = @(t) ppval(Upp,t);

% 2) Set up system dynamics
paras.v     = 1;    % constant speed
paras.r     = 3;    % truncation order for Taylor expansion
dubins      = DubinsCar(paras);
dynamics    = @dubins.dynamics;
polydyn     = @dubins.polynomialdyn;
linearize   = @dubins.linearize;
[A,B]       = linearize(x0,u0);

% 3) Compute TVLQR by solving matrix Riccati equation (backward in time!)
Q           = @(t) diag([10 10 1]);
R           = @(t) 1;
S0          = diag([2 2 4]);
tspan       = [0 t_nom(end)];
[ts,Ss]     = tvlqr_riccati(tspan,A,B,Q,R,S0);
% Spp         = spline(ts,Ss);
S           = @(t) ppval(spline(ts,Ss),t);
K           = @(t) -R(t)\B(t)'*S(t);
% sanitycheck_TVLQR;

