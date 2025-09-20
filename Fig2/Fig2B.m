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

[t,y] = ode15s(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
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

figure(4); plot(t,sT - s1bout, 'LineWidth',3); set(gca,'fontsize',35);title('Free SPEN')
ylim([0,80]); xlim([0,2000])


function dy = ODE_model(t,y,sT)
%Model parameters
a_act =  12.8; % activator synthesis rate from single X chromosome
d_act = 0.13; % degradation rate of free activator
K_n = 7.13; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 1.0; % Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 9.0; % xist synthesis rate from single X chromosome
d_x = 0.031; % degradation rate of Xist 
m = 2.0;  % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 4.3; % quantity of SPEN at which Xist dissociation is half max
K_a = 374.6; % quantity of activator at which Xist transcription is half max
k1 = 0.0020; % rate constant for Xist binding to DNA
k2 = 4.54; % maximum dissociation rate for Xist
k4 = 0.11; % dissoication rate for bound SPEN
k3 = 8.26; % association rate for SPEN
sT = sT;  % total SPEN quantity
XbsT =  100; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.

act = y(1);
x1f = y(2);
x1b = y(3);
s1b = y(4);

dy = [a_act/(1 +(s1b/K_n)^n) - d_act*act
      a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m)
      k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m)
      k3*(sT - s1b)*(N_S*x1b - s1b)  - k4*s1b
    ];
end