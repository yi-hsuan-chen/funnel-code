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
    t_nom{i}        = traj_temp(:,1);
    x_nom{i}        = traj_temp(:,2);
    y_nom{i}        = traj_temp(:,3);
    th_nom{i}       = traj_temp(:,4);
    u_nom{i}        = traj_temp(:,end);
    tGrid{i}        = traj_grid_temp(:,1);
    xGrid{i}        = traj_grid_temp(:,2);
    yGrid{i}        = traj_grid_temp(:,3);
    thGrid{i}       = traj_grid_temp(:,4);
    uGrid{i}        = traj_grid_temp(:,end);
    % Plot the entire trajectory
    plot(x_nom{i},y_nom{i},'k-','LineWidth',3);
    % Plot the grid points:
    plot(xGrid{i}, yGrid{i}, 'go','MarkerSize',3,'LineWidth',3);
    % Plot the start and end points:
    plot(x_nom{i}([1,end]), y_nom{i}([1,end]),'gs','MarkerSize',12,'LineWidth',3);
end