% Time
tspan = [0:1:4000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,8]));

% Initial conditions
act = 0; s1b = 0; xf = 0; s2b = 0; s3b = 0; 
x1b = 1; x2b = 2; x3b = 3;
y0 = [act;xf;x1b;s1b;x2b;s2b;x3b;s3b];

% Solve equations
sT = 357.786636; % Total SPEN
[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
xfout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);
x2bout = y(:,5);
s2bout = y(:,6);
x3bout = y(:,7);
s3bout = y(:,8);

% Plot results
figure; plot(t, x1bout, 'LineWidth',3);  hold on; plot(t, x2bout, 'LineWidth',3); plot(t, x3bout, 'LineWidth',3); title('Bound Xist'); set(gca,'fontsize',35)
legend('ChrX #1','ChrX #2','ChrX #3','Location','best'); ylim([0,100])

figure; plot(t, s1bout, 'LineWidth',3);  hold on; plot(t, s2bout, 'LineWidth',3); plot(t, s3bout, 'LineWidth',3); title('Bound SPEN'); set(gca,'fontsize',35)
legend('ChrX #1','ChrX #2','ChrX #3','Location','best'); 

figure; plot(t, actout, 'LineWidth',3); title('Activator'); set(gca,'fontsize',35)

figure; plot(t, sT - s1bout - s2bout - s3bout, 'LineWidth',3); title('Free SPEN'); set(gca,'fontsize',35)


function dy = ODE_model(t,y,sT)

%Model parameters
a_act =  6.107926; % activator synthesis rate from single X chromosome
d_act = 0.01225; % degradation rate of free activator
K_n = 36.396536; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 4.785343; % Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 2.338358; % xist synthesis rate from single X chromosome
d_x = 0.021799; % degradation rate of Xist 
m = 1.686465;  % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 2.621157; % quantity of SPEN at which Xist dissociation is half max
K_a = 2536.207138; % quantity of activator at which Xist transcription is half max
k1 = 0.001165; % rate constant for Xist binding to DNA
k2 = 7.452107; % maximum dissociation rate for Xist
k4 = 0.135758; % association rate for SPEN
k3 = 8.270818; % dissoication rate for bound SPEN
sT = sT;  % total SPEN quantity
XbsT = 100.513847; % quantity of Xist binding sites
N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist. 
                      % We let each chromosome be able to recruit and bind to all SPEN.
act = y(1);
xf = y(2);
x1b = y(3);
s1b = y(4);

x2b = y(5);
s2b = y(6);

x3b = y(7);
s3b = y(8);

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) + a_act/(1 +(s3b/K_n)^n) - d_act*act
      a_x*act/(K_a + act) - d_x*xf - k1*(XbsT - x1b)*xf - k1*(XbsT - x2b)*xf - k1*(XbsT - x3b)*xf + k2*x1b/(1+(s1b/K_S)^m) + k2*x2b/(1+(s2b/K_S)^m) + k2*x3b/(1+(s3b/K_S)^m)
      k1*(XbsT - x1b)*xf - k2*x1b/(1+(s1b/K_S)^m)
      k3*(sT - s1b - s2b - s3b)*(N_S*x1b - s1b)  - k4*s1b
      k1*(XbsT - x2b)*xf - k2*x2b/(1+(s2b/K_S)^m)
      k3*(sT - s1b - s2b - s3b)*(N_S*x2b - s2b) - k4*s2b
      k1*(XbsT - x3b)*xf - k2*x3b/(1+(s3b/K_S)^m)
      k3*(sT - s1b - s2b - s3b)*(N_S*x3b - s3b)  - k4*s3b
    ]; 
end