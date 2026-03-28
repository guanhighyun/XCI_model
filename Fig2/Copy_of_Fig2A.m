% Time
tspan = [0,2000];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,7]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; 
x1b = 1; x2b = 2;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b];

% Solve equations
sT = 105; % Total SPEN
[t,y] = ode15s(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
x1bout = y(:,3);
s1bout = y(:,4);
x2bout = y(:,6);
s2bout = y(:,7);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  hold on; plot(t, x2bout, 'LineWidth',3); ylim([0,100]); set(gca,'fontsize',35); title('Bound Xist')
legend('ChrX #1','ChrX #2','location','best')

figure; plot(t, s1bout, 'LineWidth',3);  hold on; plot(t, s2bout, 'LineWidth',3); set(gca,'fontsize',35); title('Bound SPEN')
legend('ChrX #1','ChrX #2','location','best')

figure; plot(t, actout, 'LineWidth',3);  set(gca,'fontsize',35); title('Activator')

figure; plot(t, sT - s1bout - s2bout, 'LineWidth',3); set(gca,'fontsize',35); title('Free SPEN')

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

x2f = y(5);
x2b = y(6);
s2b = y(7);

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) - d_act*act;
    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x1b - s1b)  - k4*s1b;
    a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)^m);
    k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)^m);
    k3*(sT - s1b - s2b)*(N_S*x2b - s2b)  - k4*s2b;
    ]; 
end