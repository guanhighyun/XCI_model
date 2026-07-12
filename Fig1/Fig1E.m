clear;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
% Time
tspan = [0:0.1:1500];

% Error tolerance
options = odeset('RelTol',1e-13,'AbsTol',repmat(1e-13,[1,10]));

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0;
x1b = 1; x2b = 2; x3b = 10;
y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b];

% Total SPEN
sT = 1000;

N_c = 3;
p = xci_params('activator_inhibition');

% Solve equations
[t,y] = ode45(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
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
legend('ChrX #1','ChrX #2','ChrX #3','Location','best'); xlim([0,1500]); ylim([0,100])

figure, plot(t, s1bout, 'LineWidth',3),  hold on, plot(t, s2bout, 'LineWidth',3), plot(t, s3bout, 'LineWidth',3), title('Bound SPEN'), set(gca,'fontsize',35)
xlim([0,1500]); ylim([0,25]); legend('ChrX #1','ChrX #2','ChrX #3','Location','best');

figure, plot(t, actout, 'LineWidth',3), title('Activator'), set(gca,'fontsize',35)
xlim([0,1500]); ylim([0,5])

figure, plot(t, sT - s1bout - s2bout - s3bout, 'LineWidth',3), title('Free SPEN'), set(gca,'fontsize',35)
ylim([0,Inf])