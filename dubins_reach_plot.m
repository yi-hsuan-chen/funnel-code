clc; clear; close all;

projDims = [1 2];
figure; hold on; grid on;

for k = 1:5
filename    = ['reachableSet_' num2str(k,'%d')];
load(filename);

% plot reachable sets
plot(R,projDims,'Order',3,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);

% plot initial set
plot(initSet,projDims,'w','FaceColor',[.8 .8 .8],'LineWidth',2);

% plot reference trajectory
plot(xyRef(1,:),xyRef(2,:),'b','LineWidth',2);

end 

xlabel('$x$'); ylabel('$y$');
set(gca().XLabel,'Interpreter','latex');
set(gca().YLabel,'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',14);