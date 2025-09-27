% Perform scan of different initial conditions of bound Xist
% for XX case

% Simulation time points
tspan = [0:50:6000];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-10,[1,7]));

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
        [t,y] = ode15s(@(tt,yy) ODE_model(tt,yy),tspan,y0,options);
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
regime(x2bout>=50) = regime(x2bout>=50) + 1

regime(x2bout<50 & x1bout<50) = regime(x2bout<50 & x1bout<50) + 2;;

% Set aside simulations with symmetric initial conditions
regime(find(eye(size(regime)))) = 5;
figure; imagesc(regime)

%xlabel('Initial Xist on X-1')
%ylabel('Initial Xist on X-2')
set(gca,'ydir','normal')
set(gca,'fontsize',25)
axis square
caxis([0,5])

function dy = ODE_model(t,y)

%Model parameters
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

% a_act =  0.054; % activator synthesis rate from single X chromosome
% d_act = 0.00846; % degradation rate of free activator 
% K_n = 1.0; % quantity of bound SPEN at which activator synthesis rate is half max.
% n = 9.9;% Hill coefficient for SPEN supressing activator synthesis rate. 
% a_x = 6.05; % xist synthesis rate from single X chromosome
% d_x = 0.0090; % degradation rate of Xist 
% m = 3.18; % Hill coefficient for bound SPEN reducing dissociation rate Xist
% K_S = 3.81; % quantity of SPEN at which Xist dissociation is half max
% K_a = 7.73; % quantity of activator at which Xist transcription is half max
% k1 = 0.0016; % rate constant for Xist binding to DNA
% k2 = 1.78; % maximum dissociation rate for Xist
% k4 = 5.5; % dissoication rate for bound SPEN
% k3 = 0.091; % association rate for SPEN
% sT = 51;  % total SPEN quantity
% XbsT =  100; % quantity of Xist binding sites
% N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
%                       % We let each chromosome be able to recruit and bind to all SPEN.

act = y(1); % Xist activator

x1f = y(2); % free Xist produced by chromosome X1
x1b = y(3); % bound Xist on chromosome X1
s1b = y(4); % bound SPEN on chromosome X1

x2f = y(5); % free Xist produced by chromosome X2
x2b = y(6); % bound Xist on chromosome X2
s2b = y(7); % bound SPEN on chromosome X2

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) - d_act*act;

    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x1b - s1b) - k4*s1b;

    a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)^m);
    k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x2b - s2b) - k4*s2b;
    ]; 
end