
% Time
tspan = [0:0.1:5000];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,4]));

% Initial conditions
act = 0; x1f = 0;
x1b = 1; s1b = 0;
y0 = [act;x1f;x1b;s1b];

% Solve equations
sT = 135; % Total SPEN

[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
x1fout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  title('Bound Xist'); 
ylim([0,100]); xlim([0,5000]); set(gca,'fontsize',35)

figure; plot(t, s1bout, 'LineWidth',3); title('Bound SPEN'); 
ylim([0,15]); xlim([0,5000]); set(gca,'fontsize',35)

figure; plot(t, actout, 'LineWidth',3); title('Activator'); 
ylim([0,80]); xlim([0,5000]); set(gca,'fontsize',35)

% figure; plot(t, sT - s1bout, 'LineWidth',3); title('Free SPEN'); xlim([0,3500])
% ylim([0,1000]); xlim([0,1500]); set(gca,'fontsize',35)

function dy = ODE_model(t,y,sT)

%Model parameters
a_act =  0.27; % activator synthesis rate from single X chromosome
d_act = 0.0068; % degradation rate of free activator 
K_n = 1.0; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 9.63;% Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 0.35; % xist synthesis rate from single X chromosome
d_x = 0.0091; % degradation rate of Xist 
m = 2.84; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 1.65; % quantity of SPEN at which Xist dissociation is half max
K_a = 1.83; % quantity of activator at which Xist transcription is half max
k1 = 0.0082; % rate constant for Xist binding to DNA
k2 = 2.33; % maximum dissociation rate for Xist
k4 = 0.28; % dissoication rate for bound SPEN
k3 = 0.00027; % association rate for SPEN
sT = sT;  % total SPEN quantity
XbsT =  100; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.

act = y(1); % free activator 
x1f = y(2); % free Xist
x1b = y(3); % bound Xist
s1b = y(4); % bound SPEN

dy = [a_act/(1 +(s1b/K_n)^n) - d_act*act;

    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k3*(sT - s1b)*(N_S*x1b - s1b) - k4*s1b;
    ]; 
end