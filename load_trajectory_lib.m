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