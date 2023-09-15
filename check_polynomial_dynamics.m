clc; clear; close all;

load_trajectory_lib; close
clearvars -except t_nom x_nom y_nom th_nom u_nom;
Upp         = spline(t_nom,u_nom);
u0          = @(t) ppval(Upp,t);

paras.v     = 1;    % constant speed
paras.r     = 3;    % truncation order for Taylor expansion
dubins      = DubinsCar(paras);

dynamics    = @dubins.dynamics;
polydyn     = @dubins.polynomialdyn;
x0          = [0;0;0];
[t,x]       = ode45(@(t,x) dynamics(t,x,u0(t)),t_nom,x0);
[~,x_app]   = ode45(@(t,x) polydyn(t,x,u0(t)),t_nom,x0);

figure();
plot(x_nom,y_nom,'k',x(:,1),x(:,2),'b--','LineWidth',2); hold on; grid on;
plot(x_app(:,1),x_app(:,2),'r-.','LineWidth',2); 
plot(x_nom([1,end]), y_nom([1,end]),'ks','MarkerSize',12,'LineWidth',3);

