addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% Time
tspan = [0:1:20000];

% Scan list of k1 and k3
k1_list = [0.0005:0.00005:0.01];
k3_list = [0.00001:0.000005:0.0007];

N_c = 4;
p = xci_params('activator_inhibition');

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,13]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;
x1b = 1; x2b = 2; x3b = 3; x4b = 4;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b;x4f;x4b;s4b];

% Solve equations
x1bout = nan(numel(k1_list),numel(k3_list));
x2bout = nan(numel(k1_list),numel(k3_list));
x3bout = nan(numel(k1_list),numel(k3_list));
x4bout = nan(numel(k1_list),numel(k3_list));

for i = 1:numel(k1_list)
    for j = 1:numel(k3_list)
        p(10) = k1_list(i);   % k1 is parameter index 10
        p(13) = k3_list(j);   % k3 is parameter index 13
        [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
        % Extract values from the 10000 to 20000 min and calculate the mean
        x1bout(i,j) = mean(y(end-10000:end,3));
        x2bout(i,j) = mean(y(end-10000:end,6));
        x3bout(i,j) = mean(y(end-10000:end,9));
        x4bout(i,j) = mean(y(end-10000:end,12));
    end
end
save('Fig3A.mat')

% Plot the bifurcation diagram
figure;
regime = zeros(numel(k1_list),numel(k3_list));

% Use a bound Xist ratio of 1.2 as the threshold to define regimes
regime(x2bout./x1bout<1.2) = 0;
regime(x2bout./x1bout>=1.2) = 3;
regime(x3bout./x1bout>=1.2 & x2bout./x1bout<2) = 2;
regime(x4bout./x1bout>=1.2 & x3bout./x1bout<2) = 1;

imagesc(regime); set(gca,'ydir','normal'); ylabel('k_1 (min^{-1})');
xlabel('k_3 (min^{-1})'); 
set(gca,'fontsize',20)
xticks([1,39:40:139])
xticklabels(k3_list([1,39:40:139]))
yticks([1,41:40:191])
yticklabels(k1_list([1,41:40:191]))