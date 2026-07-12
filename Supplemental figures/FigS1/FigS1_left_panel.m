clear;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
% Perform scan of different initial conditions of bound Xist
% for XY case

% Simulation time points
tspan = [0:10:20000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,4]));

% Initial conditions
sT = 1000;
act = 0; x1f = 0;
s1b = sT;
x1b_list = 0:100;

N_c = 1;
p = xci_params('activator_inhibition');

% Begin scan
for i = 1:numel(x1b_list)
    x1b = x1b_list(i);
    y0 = [act;x1f;x1b;s1b];

    % Solve equations
    [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
    x1bout(i) = y(end,3);
end

% Plot results
figure; plot(x1b_list,x1bout,'k.-','linewidth',4,'MarkerSize',20); ylim([0,6])
set(gca,'fontsize',25); xlabel('Initial bound Xist'); ylabel('Steady-state bound Xist')
axis square