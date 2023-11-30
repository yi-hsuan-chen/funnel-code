addpath('dynSys','lib');
clc; clear; close all;

%% Step 1: Trajectory Generation
load_trajectory_lib; close;

%% Step 2: Dynamics Setup
% 0) Set up system dynamics
paras.v     = 1;    % constant speed
paras.r     = 5;    % truncation order for Taylor expansion

dubins      = DubinsCar(paras);
linearize   = @dubins.linearize;

i = 5;
% 1) Fit (xo,uo) to piecewise polynomials
X_nom       = [x_nom{i} y_nom{i} th_nom{i}];
T_nom       = t_nom{i};
TGrid       = tGrid{i};
U_nom       = u_nom{i};

Xpp         = spline(T_nom,X_nom');
Upp         = spline(T_nom,U_nom);
x0          = @(t) ppval(Xpp,t);
u0          = @(t) ppval(Upp,t);
[A,B]       = linearize(x0,u0);

% 2) Compute TVLQR by solving matrix Riccati equation (backward in time!)
Q           = @(t) diag([10 10 1]);
R           = @(t) 1;
S0          = diag([2 2 4]);
tspan       = [0 T_nom(end)];
[ts,Ss]     = tvlqr_riccati(tspan,A,B,Q,R,S0);
S           = @(t) ppval(spline(ts,Ss),t);
K           = @(t) -R(t)\B(t)'*S(t);

%% Step 3: Compute 'input' trajectory (as a matrix)
% uMat  includes K matrix, and reference trajectory (x0,u0)
timeStep    = 0.01;
tfinal      = T_nom(end);
nSteps      = round(tfinal/timeStep);
uMat        = zeros(7,nSteps);
ts          = 0:timeStep:tfinal;
for i = 1:nSteps
    uMat(1:3,i)     = K(ts(i))';
    uMat(4:6,i)     = x0(ts(i))';
    uMat(end,i)     = u0(ts(i));
end
save('tvlqr_traj5_xend=2.mat','uMat')