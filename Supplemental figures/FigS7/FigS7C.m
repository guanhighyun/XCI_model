addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));

% Time
tspan = [0:0.1:2000];

% Error tolerance
options = odeset('RelTol',1e-10,'AbsTol',repmat(1e-10,[1,8]));

% Initial conditions
act = 0; s1b = 0; xf = 0; s2b = 0; s3b = 0; 
x1b = 1; x2b = 2; x3b = 3;
y0 = [act;xf;x1b;s1b;x2b;s2b;x3b;s3b];

% Solve equations
sT = 676; % Total SPEN
N_c = 3;
p = xci_params('figS7');
[t,y] = ode45(@(tt,yy) ODE_model_shared_xist(tt,yy,p,N_c),tspan,y0,options);
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
legend('ChrX #1','ChrX #2','ChrX #3','Location','best'); ylim([0,700])

figure; plot(t, actout, 'LineWidth',3); title('Activator'); set(gca,'fontsize',35)

figure; plot(t, sT - s1bout - s2bout - s3bout, 'LineWidth',3); title('Free SPEN'); set(gca,'fontsize',35)
