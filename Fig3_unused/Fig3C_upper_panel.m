% Time
tspan = [0:10:30000];

% Scan list of volume
V_list = [1:0.02:2];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,13]));

%initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;

% Solve equations
x1bout = nan(1,numel(V_list));
x2bout = nan(1,numel(V_list));
x3bout = nan(1,numel(V_list));
x4bout = nan(1,numel(V_list));

flag = 1;
for i = 1:numel(V_list)
    % Write initial conditions with units of concentrations

    if flag == 1
        x1b = 10/V_list(i); x2b = 20/V_list(i); x3b = 30/V_list(i); x4b = 40/V_list(i);
        y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b;x4f;x4b;s4b];
    
        [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,V_list(i)),tspan,y0,options);
    
        if abs(mean(y(end-100:end,12))-mean(y(end-200:end-100,12))) > 0.1
            flag = 0;
        end
    
        % Extract values from the 10000 to 20000 min and calculate the mean
        x1bout(i) = mean(y(end-1000:end,3))*V_list(i);
        x2bout(i) = mean(y(end-1000:end,6))*V_list(i);
        x3bout(i) = mean(y(end-1000:end,9))*V_list(i);
        x4bout(i) = mean(y(end-1000:end,12))*V_list(i);

    end
end
%save('Fig3C.mat')

% Plot the bifurcation diagram
figure;
regime = zeros(1,numel(V_list));

% Use a bound Xist ratio of 1.2 as the threshold to define regimes
regime(x2bout./x1bout<1.2) = 0;
regime(x2bout./x1bout>=1.2) = 3;
regime(x3bout./x1bout>=1.2 & x2bout./x1bout<2) = 2;
regime(x4bout./x1bout>=1.2 & x3bout./x1bout<2) = 1;

imagesc(regime); set(gca,'ydir','normal');
xlabel('V'); xticks([1,10:10:numel(V_list)]); xticklabels([1,1.2:0.2:2])
yticks([]); set(gca, 'TickDir', 'out');

function dy = ODE_model(t,y,V)

a_act =  0.27/V; % activator synthesis rate from single X chromosome
d_act = 0.0068; % degradation rate of free activator 
K_n = 1.0; % quantity of bound SPEN at which activator synthesis rate is half max.
n = 9.63;% Hill coefficient for SPEN supressing activator synthesis rate. 
a_x = 0.35/V; % xist synthesis rate from single X chromosome
d_x = 0.0091; % degradation rate of Xist 
m = 2.84; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 1.65; % quantity of SPEN at which Xist dissociation is half max
K_a = 1.83; % quantity of activator at which Xist transcription is half max
k1 = 0.0082*V; % rate constant for Xist binding to DNA
k2 = 2.33; % maximum dissociation rate for Xist
k4 = 0.28; % dissoication rate for bound SPEN
k3 = 0.00027*V; % association rate for SPEN
sT = 135/V;  % total SPEN quantity
XbsT =  100/V; % quantity of Xist binding sites
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
    k3*(sT - s1b - s2b - s3b - s4b)*(N_S*x1b - s1b)  - k4*s1b;

    a_x*act/(K_a + act) - d_x*x2f - k1*(XbsT - x2b)*x2f + k2*x2b/(1+(s2b/K_S)^m);
    k1*(XbsT - x2b)*x2f - k2*x2b/(1+(s2b/K_S)^m);
    k3*(sT - s1b - s2b - s3b - s4b)*(N_S*x2b - s2b)  - k4*s2b;

    a_x*act/(K_a + act) - d_x*x3f - k1*(XbsT - x3b)*x3f + k2*x3b/(1+(s3b/K_S)^m);
    k1*(XbsT - x3b)*x3f - k2*x3b/(1+(s3b/K_S)^m);
    k3*(sT - s1b - s2b - s3b - s4b)*(N_S*x3b - s3b)  - k4*s3b;

    a_x*act/(K_a + act) - d_x*x4f - k1*(XbsT - x4b)*x4f + k2*x4b/(1+(s4b/K_S)^m);
    k1*(XbsT - x4b)*x4f - k2*x4b/(1+(s4b/K_S)^m);
    k3*(sT - s1b - s2b - s3b - s4b)*(N_S*x4b - s4b)  - k4*s4b;
    ]; 
end