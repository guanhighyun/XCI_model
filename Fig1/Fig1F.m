close all;
rng(1)
% Load optimized model parameters from EA
param=readtable("../Parameters_for_activator_inhibition.csv");

% Plot distribution of m
swarmchart(ones(1,numel(param.m)),param.m,100,'filled'); xlim([-1,3]);
xticks([]); set(gca,'fontsize',25); ylim([0,7]); 
hold on; plot([-5,7],[1,1],'r','LineWidth',5);
ylabel('Value'); set(gca,'fontsize',35)
