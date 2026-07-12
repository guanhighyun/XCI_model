close all; clear;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
% Time
tspan = [0:0.1:200];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,4]));

% Initial conditions
act = 0; x1f = 0;
x1b = 1; s1b = 0;
y0 = [act;x1f;x1b;s1b];

% Solve equations
sT = 1303; % Total SPEN

N_c = 1;
p = xci_params('figS4');

[t,y] = ode45(@(tt,yy) ODE_model_promoted_binding(tt,yy,p,N_c,true),tspan,y0,options);
actout = y(:,1);
x1fout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  title('Bound Xist'); xlim([0,3500])
ylim([0,124]); xlim([0,200]); set(gca,'fontsize',35)

figure; plot(t, s1bout, 'LineWidth',3); title('Bound SPEN'); xlim([0,3500])
ylim([0,1400]); xlim([0,200]); set(gca,'fontsize',35)

figure; plot(t, actout, 'LineWidth',3); title('Activator'); xlim([0,3500])
ylim([0,350]); xlim([0,200]); set(gca,'fontsize',35)

figure; plot(t, sT - s1bout, 'LineWidth',3); title('Free SPEN'); xlim([0,3500])
ylim([0,1400]); xlim([0,200]); set(gca,'fontsize',35)

