addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..'));

% Perform scan of different initial conditions of bound Xist
% for XXX case

% Simulation time points
tspan = [0:10:20000];

% Initial conditions
act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0;

x1b_IC = 0:2:100; 
x2b_IC = 0:2:100; 
x3b_IC = 0:20:100; 

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-10,[1,10]));

x1bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));
s1bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));
x2bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));
s2bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));
x3bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));
s3bout = nan(numel(x1b_IC), numel(x2b_IC),numel(x3b_IC));

N_c = 3;
p = xci_params('fig2');

% Begin scan
for i = 1:numel(x1b_IC)
    for j = 1:numel(x2b_IC)
        for k = 1:numel(x3b_IC)
            y0 = [act;x1f;x1b_IC(i);s1b;x2f;x2b_IC(j);s2b;x3f;x3b_IC(k);s3b];
            [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,p,N_c),tspan,y0,options);
            x1bout(i,j,k) = y(end,3);
            s1bout(i,j,k) = y(end,4);
            x2bout(i,j,k) = y(end,6);
            s2bout(i,j,k) = y(end,7);
            x3bout(i,j,k) = y(end,9);
            s3bout(i,j,k) = y(end,10);
        end
    end
end

% Plot bifurcation diagram with regard to different initial donditions of
% bound Xist. Use 50 as the threshold value to define successful
% establishment of an inactive X.
regime = zeros(numel(x1b_IC), numel(x2b_IC), numel(x3b_IC));
regime(x1bout>=50) = regime(x1bout>=50) + 1;
regime(x2bout>=50) = regime(x2bout>=50) + 1;
regime(x3bout>=50) = regime(x3bout>=50) + 1;


% Set aside simulations with symmetric initial conditions
for i = 1:numel(x1b_IC)
    for j = 1:numel(x2b_IC)
        for k = 1:numel(x3b_IC)
            if ((x1b_IC(i) == x2b_IC(j)) || (x2b_IC(j) == x3b_IC(k)) || (x1b_IC(i) == x3b_IC(k)))
                regime(i,j,k) = 5;
            end
        end
    end
end

figure; 
h = slice(x1b_IC,x2b_IC,x3b_IC, regime,[],[],0:20:100);
set(gca,'fontsize',20)
set(h,'edgecolor','none');
xlabel({''})
ylabel({''})
zlabel({''})

xlim([0,100])
ylim([0,100])
zlim([0,100])
view(-45,45)
grid off; axis square;
alpha(0.8)
caxis([0,5])
