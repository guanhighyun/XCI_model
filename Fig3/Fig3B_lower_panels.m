clear;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
% Time
tspan = [0:0.5:600];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-8,[1,13]));

% Set the value of nuclear volumes
V_list = [1, 1.3, 1.5, 2];

N_c = 4;
p = xci_params('activator_inhibition');

for i = 1:numel(V_list)
    V = V_list(i);
    %initial conditions
    act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;
    x1b = 1/V; x2b = 2/V; x3b = 10/V; x4b = 20/V;
    y0 = [act;x1f;x1b;s1b;x2f;x2b;s2b;x3f;x3b;s3b;x4f;x4b;s4b];

    pv = xci_scale_volume(p, V);
    pv(15) = 1000/V;      % sT in concentration units

    % Solve equations
    [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,pv,N_c),tspan,y0,options);
    actout = y(:,1)*V;
    x1fout = y(:,2)*V;
    x1bout = y(:,3)*V;
    s1bout = y(:,4)*V;
    x2fout = y(:,5)*V;
    x2bout = y(:,6)*V;
    s2bout = y(:,7)*V;
    x3fout = y(:,8)*V;
    x3bout = y(:,9)*V;
    s3bout = y(:,10)*V;
    x4fout = y(:,11)*V;
    x4bout = y(:,12)*V;
    s4bout = y(:,13)*V;

    % Plot results
    figure; plot(t, x1bout, 'LineWidth',3);  hold on; plot(t, x2bout, 'LineWidth',3); plot(t, x3bout, 'LineWidth',3); plot(t, x4bout, 'LineWidth',3); title('Bound Xist'); set(gca,'fontsize',25)
    title([]); xticks([0:200:600]); %xlim([0,600]); ylim([0,100])
    set(gca,'fontsize',40)
end