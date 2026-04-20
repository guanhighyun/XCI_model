clear;  close all;
% Time
tspan = [0:100:10000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,10]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0;
x1b = 1; x2b = 2; x3b = 10;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b];

% Total SPEN 
sT = 1676;

% Solve equations
[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,sT),tspan,y0,options);
actout = y(:,1);
x1fout = y(:,2);
x1bout = y(:,3);
s1bout = y(:,4);
x2fout = y(:,5);
x2bout = y(:,6);
s2bout = y(:,7);
x3fout = y(:,8);
x3bout = y(:,9);
s3bout = y(:,10);

% Plot figures
figure; plot(t, x1bout, 'LineWidth',3);  hold on; plot(t, x2bout, 'LineWidth',3); plot(t, x3bout, 'LineWidth',3); title('Bound Xist'); set(gca,'fontsize',35)
legend('ChrX #1','ChrX #2','ChrX #3','Location','best');  ylim([0,100]) ; 

figure, plot(t, s1bout, 'LineWidth',3),  hold on, plot(t, s2bout, 'LineWidth',3), plot(t, s3bout, 'LineWidth',3), title('Bound SPEN'), set(gca,'fontsize',35)
 legend('ChrX #1','ChrX #2','ChrX #3','Location','best');

ylim([0,20])

figure, plot(t, actout, 'LineWidth',3), title('Activator'), set(gca,'fontsize',35)
ylim([0,10]); 

figure, plot(t, sT - s1bout - s2bout - s3bout, 'LineWidth',3), title('Free SPEN'), set(gca,'fontsize',35)
ylim([0,700]);

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

act = y(1); % Xist activator

x1f = y(2); % free Xist produced by chromosome X1
x1b = y(3); % bound Xist on chromosome X1
s1b = y(4); % bound SPEN on chromosome X1

x2f = y(5); % free Xist produced by chromosome X2
x2b = y(6); % bound Xist on chromosome X2
s2b = y(7); % bound SPEN on chromosome X2

x3f = y(8); % free Xist produced by chromosome X3
x3b = y(9); % bound Xist on chromosome X3
s3b = y(10); % bound SPEN on chromosome X3

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) + a_act/(1 +(s3b/K_n)^n) - d_act*act
        a_x*act/(K_a + act) - d_x*x1f - (XbsT - x1b)*x1f*(k1+k5*(s1b^m)/(s1b^m + K_S^m)) + k2*x1b
        (XbsT - x1b)*x1f*(k1+k5*(s1b^m)/(s1b^m + K_S^m)) - k2*x1b
        k3*(sT - s1b - s2b - s3b)*(N_S*x1b - s1b)  - k4*s1b
        a_x*act/(K_a + act) - d_x*x2f - (XbsT - x2b)*x2f*(k1+k5*(s2b^m)/(s2b^m + K_S^m)) + k2*x2b
        (XbsT - x2b)*x2f*(k1+k5*(s2b^m)/(s2b^m + K_S^m)) - k2*x2b
        k3*(sT - s1b - s2b - s3b)*(N_S*x2b - s2b)  - k4*s2b
        a_x*act/(K_a + act) - d_x*x3f - (XbsT - x3b)*x3f*(k1+k5*(s3b^m)/(s3b^m + K_S^m)) + k2*x3b
        (XbsT - x3b)*x3f*(k1+k5*(s3b^m)/(s3b^m + K_S^m)) - k2*x3b
        k3*(sT - s1b - s2b - s3b)*(N_S*x3b - s3b)  - k4*s3b
    ]; 
end