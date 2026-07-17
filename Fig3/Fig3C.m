clear;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
% Time
tspan = [0:1:20000];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-8,[1,13]));

N_c = 4;
p = xci_params('activator_inhibition');

% Set the value of nuclear volumes (Up to 2X volume as per manuscript)
V_list = 1:0.002:2; 
sT = 1000;
new_sT_list = sT * (1:0.002:2);

x1bout = nan(numel(V_list),numel(new_sT_list));
x2bout = nan(numel(V_list),numel(new_sT_list));
x3bout = nan(numel(V_list),numel(new_sT_list));
x4bout = nan(numel(V_list),numel(new_sT_list));

for i = 1:numel(V_list)
    V = V_list(i);
   
    
    for j = 1:numel(new_sT_list)
        new_sT = new_sT_list(j); 
        
        % Initial conditions
        act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;
        x1b = 1; x2b = 2; x3b = 3; x4b = 4;
        y0 = [act; x1f; x1b; s1b; x2f; x2b; s2b; x3f; x3b; s3b; x4f; x4b; s4b];

        pv = xci_scale_volume(p, V);
        pv(15) = new_sT/V;    % sT in concentration units

        % Solve equations
        [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,pv,N_c), tspan, y0, options);
        
        % Extract values from the 10000 to 20000 min and calculate the mean
        x1bout(i,j) = mean(y(end-10000:end,3));
        x2bout(i,j) = mean(y(end-10000:end,6));
        x3bout(i,j) = mean(y(end-10000:end,9));
        x4bout(i,j) = mean(y(end-10000:end,12));
        
   end
end
save('Fig3C.mat')

figure;
regime = zeros(numel(V_list),numel(new_sT_list));

% Use a bound Xist ratio of 1.2 as the threshold to define regimes
regime(x2bout./x1bout < 1.2) = 0;
regime(x2bout./x1bout >= 1.2) = 3;
regime(x3bout./x1bout >= 1.2 & x2bout./x1bout < 2) = 2;
regime(x4bout./x1bout >= 1.2 & x3bout./x1bout < 2) = 1;

% Plot
imagesc(V_list, new_sT_list/sT, regime');
set(gca,'YDir','normal');

xlabel('Relative volume');
ylabel('Relative SPEN abundance S/S^T');
set(gca,'FontSize',20);

xticks(1:0.2:2);
yticks(1:0.2:2);