addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));
% Perform scan of different initial conditions of bound Xist
% for XX case

% Simulation time points
tspan = [0:50:6000];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-10,[1,7]));

N_c = 2;
p = xci_params('activator_inhibition');

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0;

% List of bound Xist initial conditions
x1b_IC = 0:1:100;
x2b_IC = 0:1:100;

x1bout = nan(numel(x1b_IC),numel(x2b_IC));
s1bout = nan(numel(x1b_IC),numel(x2b_IC));
x2bout = nan(numel(x1b_IC),numel(x2b_IC));
s2bout = nan(numel(x1b_IC),numel(x2b_IC));

% Perform scan
for i = 1:numel(x1b_IC)
    for j = 1:numel(x2b_IC)
        y0 = [act;x1f;x1b_IC(i);s1b;x2f;x2b_IC(j);s2b];
        [t,y] = ode15s(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
        x1bout(i,j) = y(end,3);
        s1bout(i,j) = y(end,4);
        x2bout(i,j) = y(end,6);
        s2bout(i,j) = y(end,7);
    end
end

% Plot bifurcation diagram with regard to different initial donditions of
% bound Xist. Use 50 as the threshold value to define successful
% establishment of an inactive X.
regime = zeros(numel(x1b_IC), numel(x2b_IC));
regime(x1bout>=50) = regime(x1bout>=50) + 1;
regime(x2bout>=50) = regime(x2bout>=50) + 1;

% Set aside simulations with symmetric initial conditions
regime(find(eye(size(regime)))) = 5;
imagesc(regime)

xlabel('Initial Xist on X-1')
ylabel('Initial Xist on X-2')
set(gca,'ydir','normal')
set(gca,'fontsize',25)
axis square
caxis([0,5])