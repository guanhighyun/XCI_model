close all; clear;
% Time
tspan = [0:0.1:800];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,4]));

% Initial conditions
act = 0; x1f = 0;
x1b = 1; s1b = 0;
y0 = [act;x1f;x1b;s1b];

% Solve equations
sT = 676; % Total SPEN

[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
x1fout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  title('Bound Xist'); ylim([0,100]) ; xlim([0,800]); set(gca,'fontsize',35)

figure; plot(t, s1bout, 'LineWidth',3); title('Bound SPEN'); xlim([0,3500])
 set(gca,'fontsize',35);  xlim([0,800])
ylim([0,20])

figure; plot(t, actout, 'LineWidth',3); title('Activator'); xlim([0,3500])
ylim([0,10]); xlim([0,800]); set(gca,'fontsize',35)

figure; plot(t, sT - s1bout, 'LineWidth',3); title('Free SPEN'); xlim([0,3500])
ylim([0,700]); xlim([0,800]); set(gca,'fontsize',35)

function dy = ODE_model(t,y,sT)

%Model parameters
a_act =  1.89; % activator synthesis rate from single X chromosome
d_act = 3.43; % degradation rate of free activator 
K_n = 2.85; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 7.8;% Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 10.0; % xist synthesis rate from single X chromosome
d_x = 0.029; % degradation rate of Xist 
m = 3.6; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 5.78; % quantity of SPEN at which Xist-chromosome binding rate is half max
K_a = 4.15; % quantity of activator at which Xist transcription is half max
k1 = 0.011; % basal rate constant for Xist binding to DNA
k2 = 8.57; % maximum dissociation rate for Xist
k3 = 0.00019; % association rate for SPEN
k4 = 8.39; % dissoication rate for bound SPEN
k5 = 1.28; % rate at which SPEN promotes Xist-chromosome binding
sT = sT;  % total SPEN quantity
XbsT =  101; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.


act = y(1); % free activator 
x1f = y(2); % free Xist
x1b = y(3); % bound Xist
s1b = y(4); % bound SPEN

dy = [a_act/(1 +(s1b/K_n)^n) - d_act*act;

    a_x*act/(K_a + act) - d_x*x1f - (XbsT - x1b)*x1f*(k1+k5*(s1b^m)/(s1b^m + K_S^m)) + k2*x1b;
    (XbsT - x1b)*x1f*(k1+k5*(s1b^m)/(s1b^m + K_S^m)) - k2*x1b;
    k3*(sT - s1b)*(N_S*x1b - s1b) - k4*s1b;
    ]; 
end