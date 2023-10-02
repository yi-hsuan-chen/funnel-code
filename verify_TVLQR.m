fcl         = @(t,x) dynamics(t,x,u0(t)+K(t)*(x-x0(t)));
checkTVLQRctrlperf(fcl,tspan,X_nom);
[~,x_app]   = ode45(@(t,x) polydyn(t,x,u0(t)+K(t)*(x-x0(t))),tspan,[0;0;0]);
plot(x_app(:,1),x_app(:,2),'r:','LineWidth',3);

function checkTVLQRctrlperf(fcl,tspan,X_nom)
% Note: Check control performance of TVLQR for nonlinear sys
x0          = [0;0;0];
[~,x_n]     = ode45(fcl,tspan,x0);
x_nom       = X_nom(:,1);
y_nom       = X_nom(:,2);
plot(x_n(:,1),x_n(:,2),'b--','LineWidth',3); grid on; 
end