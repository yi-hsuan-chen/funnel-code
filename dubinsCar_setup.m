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
    verify_TVLQR;

    %% Step 3: Implement Algorithm 1 (This should be in the while loop)

    % 3.1) Initialization: V0,V0dot,rho,rhodot,Phi,ui
    % ui is the nominal control expressed in msspoly
    [Vtraj0,rho0pp,ui,intvs,t,x,options] = InitializeOptim(TGrid,x0,u0,K,S);

    % 3.2) Perform bilinear alternation search
    V               = Vtraj0.V0;
    rhopp           = rho0pp;
    k               = 1;
    
%     searchLangrange;
    [L1,Le,Sks]     = findLs(V,rhopp,TGrid,polydyn,options,ui,intvs,t,x,k);

    searchVRho;

        
end



for i = 1:length(Sks)
    %     det(Sks{i})
    C   = [xGrid{1}(i+1) yGrid{1}(i+1) thGrid{1}(i+1)];
    plotEllipse(Sks{i},C); hold on;
end





%%
% 




% function plotEllipse(M,C)
% N   = 50;
% P   = [1 0 0 ; 0 1 0]; % xy
% E   = inv(P*inv(M)*P');
% 
% % plots an ellipse of the form xEx = 1
% R   = chol(E^(1/2));
% t   = linspace(0, 2*pi, N); % or any high number to make curve smooth
% z   = [cos(t); sin(t)];
% ell = inv(R)*z;
% 
% fill(C(1) + ell(1,:),C(2) + ell(2,:),[0.8 0.8 0.8],'FaceAlpha',0.3);
% end