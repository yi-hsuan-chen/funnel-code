clc; clear; close all;

%% Step 1: Trajectory Generation
load_trajectory_lib; close;

%% Step 2: Dynamics Setup
% 0) Set up system dynamics
paras.v     = 1;    % constant speed
paras.r     = 5;    % truncation order for Taylor expansion
dubins      = DubinsCar(paras);
dynamics    = @dubins.dynamics;
polydyn     = @dubins.polynomialdyn;
linearize   = @dubins.linearize;

nSim        = 5;
stepSize    = 5e-2;                                                    
testPerturbation = 0.1;  % Initial position error amplitude

for i = 1:1
    % 1) Fit (xo,uo) to piecewise polynomials
    X_nom       = [x_nom{i} y_nom{i} th_nom{i}];
    T_nom       = t_nom{i};
    TGrid       = tGrid{i};
    U_nom       = u_nom{i};

    % Plot the entire trajectory
    plot(x_nom{i},y_nom{i},'k-','LineWidth',3); hold on; grid on;
    % Plot the grid points:
    plot(xGrid{i}, yGrid{i}, 'go','MarkerSize',3,'LineWidth',3);
    % Plot the start and end points:
    plot(x_nom{i}([1,end]), y_nom{i}([1,end]),'gs','MarkerSize',12,'LineWidth',3);
%     hold off;

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
    
    for i = 1:nSim
        X0      = testPerturbation*randn(3,1);
        xbar    = @(t,x) x-x0(t);
        ubar    = @(t,x) K(t)*xbar(t,x);
        u       = @(t,x) u0(t)+ubar(t,x);
        xdot    = @(t,x) polydyn(t,x,u(t,x));
        x0dot   = @(t,x) polydyn(t,x0(t),u0(t));
        xbardot = @(t,x) polydyn(t,x0(t)+x,u0(t)+K(t)*x)-polydyn(t,x0(t),u0(t));
        [t,X]   = rk4(xdot,tspan,X0,stepSize);
        [~,Xnom]= rk4(x0dot,tspan,[0;0;0],stepSize);
        [~,Xbar]= rk4(xbardot,tspan,X0,stepSize);
%         [~,Xapp]= rk4(@(t,x) A(t)*x+B(t)*K(t)*x,tspan,X0,5e-2);
%         plot(Xapp(:,1),Xapp(:,2),'r--','LineWidth',3);
        plot(X(:,1),X(:,2),'b-','LineWidth',3);
        plot(Xnom(:,1),Xnom(:,2),'c:','LineWidth',3);
        plot(Xbar(:,1)+Xnom(:,1),Xbar(:,2)+Xnom(:,2),'r--','LineWidth',3);
%         plot(Xbar(:,1),Xbar(:,2),'LineWidth',3); hold on;
    end
%     c           = 2;
%     rhot        = exp(c*(T_nom-max(T_nom))/(max(T_nom)-min(T_nom)));
%     rho0pp      = interp1(T_nom,rhot,'linear','pp');
%     for j = 1:5:length(T_nom)
%         C   = X_nom(j,:);
%         plotEllipse(S(T_nom(j))./ppval(rho0pp,T_nom(j)),C); hold on;
%     end
end