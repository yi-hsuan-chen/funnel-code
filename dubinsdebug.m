clc; clear; close all;

%% Step 1: Trajectory Generation
load_trajectory_lib;

%% Step 2: Dynamics Setup
% 0) Set up system dynamics
paras.v     = 1;    % constant speed
paras.r     = 3;    % truncation order for Taylor expansion
dubins      = DubinsCar(paras);
dynamics    = @dubins.dynamics;
polydyn     = @dubins.polynomialdyn;
linearize   = @dubins.linearize;

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
    S0          = diag([2 2 4]);
    tspan       = [0 T_nom(end)];
    [ts,Ss]     = tvlqr_riccati(tspan,A,B,Q,R,S0);
    % Spp         = spline(ts,Ss);
    S           = @(t) ppval(spline(ts,Ss),t);
    K           = @(t) -R(t)\B(t)'*S(t);
%     verify_TVLQR;

    %% Step 3: Implement Algorithm 1

    % 3.1) Initialization: V0,V0dot,rho,rhodot,Phi,ui
    % ui is the nominal control expressed in msspoly
    [Vtraj0,rho0pp,ui,intvs,t,x,options] = InitializeOptim(TGrid,x0,u0,K,S);


    % 3.2) Perform bilinear alternation search
    V               = Vtraj0.V0;
    rhopp           = rho0pp;
    k               = 1;
    
    [L1,Le]         = findLs(V,rhopp,TGrid,polydyn,options,ui,intvs,t,x,k);

    % Search for V and rho
        N           = length(TGrid);
        disp('Step 2: Searching for V and rho...');
        % Initialize SOS program
        prog        = spotsosprog;
        prog        = prog.withIndeterminate([t;x]);
        for i = 1:N-1
            [prog,rhos{i}]  = prog.newFreePoly(monomials(t,0:1));
        end
        dts         = diff(ts);
        % Build rhos and rhodots
        rhoss       = repmat(rhos,1,k);
        rhodotss    = repmat(rhos,1,k);
        for i = 1:N-1
            tt          = linspace(intvs{i}(1),intvs{i}(2),k);
            rhoss{i}    = msubs(rhos{i},t,tt);
            rhodotss{i} = msubs(diff(rhos{i},t),t,tt);
        end
        [prog,Vs{i}]    = prog.newFreePoly(monomials(t,0:1));
end
