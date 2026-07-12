addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% Time
tspan = [0:1:20000];

% Scan list of volume
V_list = [1:0.02:2];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,13]));

N_c = 4;
p = xci_params('activator_inhibition');

%initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;

% Solve equations
x1bout = nan(1,numel(V_list));
x2bout = nan(1,numel(V_list));
x3bout = nan(1,numel(V_list));
x4bout = nan(1,numel(V_list));

for i = 1:numel(V_list)
    V = V_list(i);

    % Write initial conditions with units of concentrations
    %x1b = 1/V_list(i); x2b = 2/V_list(i); x3b = 3/V_list(i); x4b = 4/V_list(i);
    x1b = 1; x2b = 2; x3b = 3; x4b = 4;
    y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b;x4f;x4b;s4b];

    pv = xci_scale_volume(p, V);
    pv(15) = 1000/V;      % sT in concentration units

    [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,pv,N_c),tspan,y0,options);

    % Extract values from the 10000 to 20000 min and calculate the mean
    x1bout(i) = mean(y(end-10000:end,3))*V_list(i);
    x2bout(i) = mean(y(end-10000:end,6))*V_list(i);
    x3bout(i) = mean(y(end-10000:end,9))*V_list(i);
    x4bout(i) = mean(y(end-10000:end,12))*V_list(i);
end
save('Fig3C.mat')

% Plot the bifurcation diagram
figure;
regime = zeros(1,numel(V_list));

% Use a bound Xist ratio of 1.2 as the threshold to define regimes
regime(x2bout./x1bout<1.2) = 0;
regime(x2bout./x1bout>=1.2) = 3;
regime(x3bout./x1bout>=1.2 & x2bout./x1bout<2) = 2;
regime(x4bout./x1bout>=1.2 & x3bout./x1bout<2) = 1;

imagesc(regime); set(gca,'ydir','normal');
xlabel('V'); xticks([1,10:10:numel(V_list)]); xticklabels([1,1.2:0.2:2])
yticks([]); set(gca, 'TickDir', 'out');