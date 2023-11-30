clc; clear; close all;

% Reachability Settings ---------------------------------------------------
options.timeStep        = 0.01;
options.taylorTerms     = 5;
options.zonotopeOrder   = 200;
options.alg             = 'lin';
options.tensorOrder     = 2;

% Parameters ---------------------------------------------------
% intial set
initSet         = zonotope([[0; 0; 0],0.15*diag([1, 1, 1])]);
params.R0       = initSet;

% reference trajectory
load('tvlqr_traj5_xend=2.mat');
params.u        = uMat;
xyRef           = [params.u(4,:);params.u(5,:)];

params.tFinal   = size(uMat,2)*options.timeStep;


% uncertain inputs
% params.U        = zonotope([0*params.u(:,1), 0.05*eye(7)]);


% System Dynamics ---------------------------------------------------------
vehicle = nonlinearSys(@dubins_model_controlled,3,7);

% Reachability Analysis ---------------------------------------------------
tic
R = reach(vehicle, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

% Simulation settings --------------------------------------------------------------
simOpt.points = 20;
simRes = simulateRandom(vehicle, params, simOpt);

%% Visualization -----------------------------------------------------------
dims    = [1 2];
ref     = [4 5];

figure; hold on; box on
projDims = dims; projRef = ref;

% plot reachable sets
plot(R,projDims,'Order',3,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);

% plot initial set
plot(params.R0,projDims,'w','FaceColor',[.8 .8 .8],'LineWidth',2);

% plot simulation results
plot(simRes,projDims);

% plot reference trajectory
plot(xyRef(1,:),xyRef(2,:),'b','LineWidth',2);

% label plot
xlabel(['x_{',num2str(projDims(1)),'}']);
ylabel(['x_{',num2str(projDims(2)),'}']);

save('reachableSet_5.mat','R','xyRef','initSet');