addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% Time
tspan = [0,2000];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-10,[1,4]));

% Initial conditions
act = 0; s1b = 0; x1f = 0;
x1b = 1;
y0 = [act;x1f;x1b;s1b];

% Solve equations

% Total SPEN 
sT = 70;

N_c = 1;
p = xci_params('fig2');

[t,y] = ode15s(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
actout = y(:,1);
x1bout = y(:,3);
s1bout = y(:,4);

% Plot results
figure(1); plot(t, x1bout, 'LineWidth',3); hold on; title("Bound Xist")
set(gca,'fontsize',35); ylim([0,100]); xlim([0,2000])

figure(2); plot(t, s1bout, 'LineWidth',3); hold on;  title("Bound SPEN")
set(gca,'fontsize',35); ylim([0,80]); xlim([0,2000])

figure(3); plot(t,actout, 'LineWidth',3); set(gca,'fontsize',35);title('Activator')
ylim([0,200]); xlim([0,2000])

% figure(4); plot(t,sT - s1bout, 'LineWidth',3); set(gca,'fontsize',35);title('Free SPEN')
% ylim([0,80]); xlim([0,2000])