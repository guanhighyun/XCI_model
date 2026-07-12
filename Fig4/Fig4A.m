addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% Time
tspan = [0:0.1:9000];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,4]));

% Initial conditions
act = 0; x1f = 0;
x1b = 1; s1b = 0;
y0 = [act;x1f;x1b;s1b];

% Solve equations
sT = 250; % Total SPEN

N_c = 1;
p = xci_params('spen_depletion');

[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
actout = y(:,1);
x1fout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  title('Bound Xist'); 
ylim([0,100]); xlim([0,9000]); set(gca,'fontsize',35)

figure; plot(t, s1bout, 'LineWidth',3); title('Bound SPEN'); 
ylim([0,80]); xlim([0,9000]); set(gca,'fontsize',35)

figure; plot(t, actout, 'LineWidth',3); title('Activator'); 
ylim([0,1300]); xlim([0,9000]); set(gca,'fontsize',35)

% figure; plot(t, sT - s1bout, 'LineWidth',3); title('Free SPEN'); xlim([0,3500])
% ylim([0,1000]); xlim([0,2000]); set(gca,'fontsize',35)
