clear;
rng(1)
load('FigS5B_input.mat');

figure; hold on;

% Time span
dt = 0.001; 
tspan = 0:dt:6000; 
nSteps = length(tspan);

% Noise intensity
sigma = 1; 

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0;

% Number of simulations
N_sim = 10;
colors = [0.4,0.8,1; 0.2,0,1; 1,0,1];

% Storage
final_x1b = cell(numel(row), N_sim);
final_x2b = cell(numel(row), N_sim);

idx_a = cell(numel(row),1);
idx_b = cell(numel(row),1);
idx_c = cell(numel(row),1);

TOTAL_A = 0; TOTAL_B = 0; TOTAL_C = 0;

for k = 1:numel(row)

    x1b = x1b_IC(row(k)); 
    x2b = x2b_IC(col(k)); 
    y0 = [act; x1f; x1b; s1b; x2f; x2b; s2b];

    for j = 1:N_sim
        Y = zeros(length(y0), nSteps);
        Y(:,1) = y0;

        traj_x1b = nan(1,6000);
        traj_x2b = nan(1,6000);

        for i = 2:nSteps
            y = Y(:,i-1);

            % Stochastic part
            g = sigma * eye(length(y0));

            % Deterministic part
            f = ODE_model(y);

            % Euler-Maruyama step
            dW = sqrt(dt) * randn(length(y0),1);
            new_Y = y + f * dt + g * dW;
            new_Y(new_Y<0) = 0;
            Y(:,i) = new_Y;

            % Save only at integer times
            if mod(i*dt,1) == 0
                idx = int32(i*dt);
                traj_x1b(idx) = new_Y(3);
                traj_x2b(idx) = new_Y(6);
            end
        end

        % Store whole trajectories at integer times
        final_x1b{k,j} = traj_x1b;
        final_x2b{k,j} = traj_x2b;

        % Classify outcomes
        if traj_x1b(end) >= 50
            idx_a{k} = [idx_a{k}, j];
            TOTAL_A = TOTAL_A + 1;
        elseif traj_x2b(end) >= 50
            idx_b{k} = [idx_b{k}, j];
            TOTAL_B = TOTAL_B + 1;
        else
            idx_c{k} = [idx_c{k}, j];
            TOTAL_C = TOTAL_C + 1;
        end
    end

    % ---- Plot results ----
    t = 1:6000;

   
end


for k = 1:numel(row)
 % Group A
    for j = idx_a{k}
        plot(t, final_x1b{k,j}, 'Color', colors(1,:), 'LineWidth', 0.5);
    end

    % Group B
    for j = idx_b{k}
        plot(t, final_x2b{k,j}, 'Color', colors(1,:), 'LineWidth', 0.5);
    end

    % Group C (neither)
     for j = idx_c{k}
        plot(t, final_x1b{k,j}, 'Color', colors(3,:), 'LineWidth', 0.5);
    end;
end


    set(gca,'FontSize',25);
    xlabel('');
    ylabel('');


function [f] = ODE_model(y)

a_act =  66.8; % activator synthesis rate from single X chromosome
d_act = 0.048; % degradation rate of free activator 
K_n = 1.0; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 10.0;% Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 5.9; % xist synthesis rate from single X chromosome
d_x = 0.0041; % degradation rate of Xist 
m = 3.85; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 7.40; % quantity of SPEN at which Xist dissociation is half max
K_a = 76.5; % quantity of activator at which Xist transcription is half max
k1 = 0.0024; % rate constant for Xist binding to DNA
k2 = 4.21; % maximum dissociation rate for Xist
k4 = 8.78; % dissoication rate for bound SPEN
k3 = 0.014; % association rate for SPEN
sT = 194;  % total SPEN quantity
XbsT =  100; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.

act = y(1);
x1f = y(2);
x1b = y(3);
s1b = y(4);

x2f = y(5);
x2b = y(6);
s2b = y(7);

f = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) - d_act*act;

    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x1b - s1b) - k4*s1b;

    a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)^m);
    k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x2b - s2b) - k4*s2b;
    ]; 
end