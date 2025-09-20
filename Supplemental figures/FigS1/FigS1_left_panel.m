clear;
% Perform scan of different initial conditions of bound Xist
% for XY case

% Simulation time points
tspan = [0:10:20000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,4]));

% Initial conditions
sT = 194;
act = 0; x1f = 0;
s1b = sT;
x1b_list = 0:100;

% Begin scan
for i = 1:numel(x1b_list)
    x1b = x1b_list(i); 
    y0 = [act;x1f;x1b;s1b];

    % Solve equations
    [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
    x1bout(i) = y(end,3);
end

% Plot results
figure; plot(x1b_list,x1bout,'k.-','linewidth',4,'MarkerSize',20); ylim([0,6])
set(gca,'fontsize',25); xlabel('Initial bound Xist'); ylabel('Steady-state bound Xist')
axis square

function dy = ODE_model(t,y,sT)

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
sT = sT;  % total SPEN quantity
XbsT =  100; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.
act = y(1);
x1f = y(2);
x1b = y(3);
s1b = y(4);

dy = [a_act/(1 +(s1b/K_n)^n) - d_act*act;

    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k3*(sT - s1b)*(N_S*x1b - s1b) - k4*s1b;
    ]; 
end