% Time
tspan = [0:0.5:2000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,6]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0;
x1b = 1; x2b = 2;
y0 = [act;x1f;x1b;s1b;x2b;s2b];

% Solve equations
sT = 676; % Total SPEN
[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
x1bout = y(:,3);
s1bout = y(:,4);
x2bout = y(:,5);
s2bout = y(:,6);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  hold on; plot(t, x2bout, 'LineWidth',3); ylim([0,100]); set(gca,'fontsize',25); title('Bound Xist')
legend('ChrX #1','ChrX #2','Location','best'); set(gca,'fontsize',35)

figure; plot(t, s1bout, 'LineWidth',3);  hold on; plot(t, s2bout, 'LineWidth',3); set(gca,'fontsize',35); title('Bound SPEN')
legend('ChrX #1','ChrX #2','Location','best'); ylim([0,400]);

figure; plot(t, actout, 'LineWidth',3);  set(gca,'fontsize',35); title('Activator'); ylim([0,80])

figure; plot(t, sT - s1bout - s2bout, 'LineWidth',3); set(gca,'fontsize',35); title('Free SPEN')
ylim([0,800])

function dy = ODE_model(t,y,sT)
%Model parameters
a_act =  0.707; % activator synthesis rate from single X chromosome
d_act = 0.028; % degradation rate of free activator
K_n = 57.8; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 4.39; % Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 3.9; % xist synthesis rate from single X chromosome
d_x = 0.027; % degradation rate of Xist 
m = 2.58;  % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 26.8; % quantity of SPEN at which Xist dissociation is half max
K_a = 351; % quantity of activator at which Xist transcription is half max
k1 = 0.0018; % rate constant for Xist binding to DNA
k2 = 1.42; % maximum dissociation rate for Xist
k4 = 0.10; % association rate for SPEN
k3 = 9.9; % dissoication rate for bound SPEN
sT = sT;  % total SPEN quantity
XbsT = 101; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.
act = y(1);
xf = y(2);
x1b = y(3);
s1b = y(4);

x2b = y(5);
s2b = y(6);

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) - d_act*act
      a_x*act/(K_a + act) - d_x*xf - k1*(XbsT - x1b)*xf - k1*(XbsT - x2b)*xf + k2*x1b/(1+(s1b/K_S)^m) + k2*x2b/(1+(s2b/K_S)^m)
      k1*(XbsT - x1b)*xf - k2*x1b/(1+(s1b/K_S)^m)
      k3*(sT - s1b - s2b)*(N_S*x1b - s1b)  - k4*s1b
      k1*(XbsT - x2b)*xf - k2*x2b/(1+(s2b/K_S)^m)
      k3*(sT - s1b - s2b)*(N_S*x2b - s2b) - k4*s2b
    ]; 
end