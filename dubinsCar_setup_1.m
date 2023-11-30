clc; clear; close all;

%% Step 1: Trajectory Generation
load_trajectory_lib;

%% Step 2: Dynamics Setup
% 0) Set up system dynamics
paras.v     = 1;    % constant speed
paras.r     = 4;    % truncation order for Taylor expansion
dubins      = DubinsCar(paras);
dynamics    = @dubins.dynamics;
polydyn     = @dubins.polydyn;
linearize   = @dubins.linearize;
% fcl         = @dubins.taylorExpandCL;

for i = 1:1
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
    Sf          = diag([2 2 4]);
    tspan       = [0 T_nom(end)];
    [ts,Ss]     = tvlqr_riccati(tspan,A,B,Q,R,Sf);
    S           = @(t) ppval(spline(ts,Ss),t);
    K           = @(t) -R(t)\B(t)'*S(t);
%     verify_TVLQR;

    %% Step 3: Implement Algorithm 1 (This should be in the while loop)
    
    % 3.1) Initialization: V0,V0dot,rho,rhodot,Phi,ui
    k               = 1;
    [Vtraj0,rho0,xdots,Phi,intvs,t,x,options] = init_Optim(TGrid,x0,u0,K,S,polydyn,k);

    % 3.2) Perform bilinear alternation search
    V               = Vtraj0.V0;
    rho             = rho0;

    converged       = 0;
    for iter = 1:3
        fprintf('=======================\n');
        fprintf('Iteration: %d \n',iter);
        
        [L1,gamma]      = findL1(V,rho,TGrid,xdots,Vtraj0,Phi,options,x);
        [VV,Rho,Phi_]   = findVRho(L1,rho,gamma,TGrid,xdots,Vtraj0,Phi,options,x);
        
        %% Reset for searchforLs
        V = VV; rho = Rho; Phi = Phi_;

    end

end

%%
for i = 1:length(Phi)

%     Sorig   = double(subs(S0s{i},t,linspace(intvs{i}(1),intvs{i}(2),1)));
    S0      = double(Vtraj0.S0s{i});
    Sk{i}   = S0+Phi{i};
    detS0   = det(S0)
    detSk   = det(Sk{i})
    C   = [xGrid{1}(i+1) yGrid{1}(i+1) thGrid{1}(i+1)];
    plotEllipse(S0,C); hold on;
    plotEllipse(Sk{i},C);
end