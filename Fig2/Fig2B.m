addpath(fullfile(fileparts(mfilename('fullpath')), '..'));

% Time
tspan = [0,2000];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,7]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; 
x1b = 1; x2b = 2;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b];

% Solve equations
sT = 70; % Total SPEN
N_c = 2;
p = xci_params('fig2');
[t,y] = ode15s(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
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
