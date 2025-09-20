% Time
tspan = [0:1:20000];

% Scan list of activator and Xist synthesis rate
a_act_list = [0.01:0.01:2];
a_x_list = [1:0.02:5];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,13]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;
x1b = 1; x2b = 2; x3b = 3; x4b = 4;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b;x4f;x4b;s4b];

% Solve equations
x1bout = nan(numel(a_act_list),numel(a_x_list));
x2bout = nan(numel(a_act_list),numel(a_x_list));
x3bout = nan(numel(a_act_list),numel(a_x_list));
x4bout = nan(numel(a_act_list),numel(a_x_list));

for i = 1:numel(a_act_list)
    for j = 1:numel(a_x_list)
        [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,a_act_list(i),a_x_list(j)),tspan,y0,options);
        % Extract values from the 10000 to 20000 min and calculate the mean
        x1bout(i,j) = mean(y(end-10000:end,3));
        x2bout(i,j) = mean(y(end-10000:end,6));
        x3bout(i,j) = mean(y(end-10000:end,9));
        x4bout(i,j) = mean(y(end-10000:end,12));
    end
end
save('Fig3A.mat')

% Plot the bifurcation diagram
figure;
regime = zeros(numel(a_act_list),numel(a_x_list));

% Use a bound Xist ratio of 1.2 as the threshold to define regimes
regime(x2bout./x1bout<1.2) = 0;
regime(x2bout./x1bout>=1.2) = 3;
regime(x3bout./x1bout>=1.2 & x2bout./x1bout<2) = 2;
regime(x4bout./x1bout>=1.2 & x3bout./x1bout<2) = 1;

imagesc(regime); set(gca,'ydir','normal'); ylabel('\alpha_{A} (min^{-1})');
xlabel('\alpha_x (min^{-1})'); 
set(gca,'fontsize',20);
xticks([1:50:201])
xticklabels(a_x_list([1:50:201]))
yticks([1,50:50:201])
yticklabels(a_act_list([1,50:50:201]))
axis square

function dy = ODE_model(t,y,a_act,a_x)
%Model parameters
a_act =  a_act; % activator synthesis rate from single X chromosome
d_act = 0.72; % degradation rate of free activator

K_n = 2.2; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 9.4; % Hill coefficient for SPEN supressing activator synthesis rate. 

a_x = a_x; % xist synthesis rate from single X chromosome
d_x = 0.0076; % degradation rate of Xist 

m = 3.09; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 2.06; % quantity of SPEN at which Xist dissociation is half max

K_a = 16.1; % quantity of activator at which Xist transcription is half max

k1 = 0.0084; % rate constant for Xist binding to DNA
k2 = 8.7; % maximum dissociation rate for Xist
k4 = 10.2; % dissoication rate for bound SPEN
k5 = 0.00026; % association rate for SPEN

sT = 1000; % total SPEN quantity
XbsT =  100; % quantity of Xist binding sites

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

x4f = y(11); % free Xist produced by chromosome X4
x4b = y(12); % bound Xist on chromosome X4
s4b = y(13); % bound SPEN on chromosome X4

dy = [a_act/(1 +(s1b/K_n)^n) + a_act/(1 +(s2b/K_n)^n) + a_act/(1 +(s3b/K_n)^n) + a_act/(1 +(s4b/K_n)^n) - d_act*act;
    
    a_x*act/(K_a + act) - d_x*x1f - k1*(XbsT - x1b)*x1f + k2*x1b/(1+(s1b/K_S)^m);
    k1*(XbsT - x1b)*x1f - k2*x1b/(1+(s1b/K_S)^m);
    k5*(sT - s1b - s2b - s3b - s4b)*(N_S*x1b - s1b)  - k4*s1b;

    a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)^m);
    k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)^m);
    k5*(sT - s1b - s2b - s3b - s4b)*(N_S*x2b - s2b)  - k4*s2b;

    a_x*act/(K_a + act) - d_x*x3f - k1*(XbsT - x3b)*x3f + k2*x3b/(1+(s3b/K_S)^m);
    k1*(XbsT - x3b)*x3f - k2*x3b/(1+(s3b/K_S)^m);
    k5*(sT - s1b - s2b - s3b - s4b)*(N_S*x3b - s3b)  - k4*s3b;

    a_x*act/(K_a + act) - d_x*x4f - k1*(XbsT - x4b)*x4f + k2*x4b/(1+(s4b/K_S)^m);
    k1*(XbsT - x4b)*x4f - k2*x4b/(1+(s4b/K_S)^m);
    k5*(sT - s1b - s2b - s3b - s4b)*(N_S*x4b - s4b)  - k4*s4b;
    ]; 
end