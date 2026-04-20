clear; 
% Time
tspan = [0:100:6000];

% Error tolerance
options = odeset('RelTol',1e-8,'AbsTol',repmat(1e-8,[1,13]));

% Set the value of nuclear volumes (Up to 2X volume as per manuscript)
V_list = 1:0.1:3; 
sT = 750;
new_sT_list = sT * (1:0.01:3.5);

% Preallocate array to store the max SPEN for each volume
spen_max = zeros(1, numel(V_list)); 

for i = 1:numel(V_list)
    V = V_list(i);
    max_tolerated_sT = NaN; % Default baseline
    found_2_of_4 = false;   % Flag to track when we enter the target regime
    
    for j = 1:numel(new_sT_list)
        new_sT = new_sT_list(j); 
        
        % Initial conditions
        act = 0; s1b = 0; x1f = 0; s2b = 0; x2f = 0; s3b = 0; x3f = 0; s4b = 0; x4f = 0;
        x1b = 1; x2b = 2; x3b = 3; x4b = 4;
        y0 = [act; x1f; x1b; s1b; x2f; x2b; s2b; x3f; x3b; s3b; x4f; x4b; s4b];
        
        % Solve equations
        [t,y] = ode45(@(tt,yy) ODE_model(tt,yy,V,new_sT), tspan, y0, options);
        
        % Extract final bound Xist values for all 4 chromosomes
        xb_final = [y(end, 3), y(end, 6), y(end, 9), y(end, 12)];
        
        % Sort descending: top 2 are the most silenced, bottom 2 are the most active
        sorted_xb = sort(xb_final, 'descend');
       
        % Check if exactly 2-of-4 silencing is maintained.
        % This is true if there is a distinct gap between the 2nd highest 
        % and 3rd highest bound Xist values.
        is_2_of_4_state = (sorted_xb(2) > sorted_xb(3) + 5);
        
        if is_2_of_4_state
            % We are in the 2-of-4 regime. Update the max tolerated value.
            max_tolerated_sT = new_sT;
            found_2_of_4 = true; % Mark that we've successfully hit the target state
            
        elseif found_2_of_4
            % If we are NO LONGER in the 2-of-4 state, but we WERE previously, 
            % it means SPEN has become too high (e.g., 3-of-4 are now silenced).
            % We have found our boundary. Break the loop to save time.
            break; 
        end
    end
    
    % Store the highest tolerated SPEN abundance for the current volume
    spen_max(i) = max_tolerated_sT;
end

% Plotting the results
figure;
plot(V_list, spen_max / sT, '-o', 'LineWidth', 4);
xlabel('Relative Nuclear Volume (V)');
ylabel('Maximum SPEN Increase (SPEN_{max} / S_T)');
title('Maximum Tolerated SPEN Abundance in Tetraploid Cells');
set(gca,'fontsize',15)
grid on; hold on; xlim([1,3.5]); ylim([1,3.5])
plot(1:0.1:3.5,1:0.1:3.5,'k--','LineWidth',3)

function dy = ODE_model(t,y,V, new_sT)
%Model parameters
a_act =  1.1/V; % activator synthesis rate from single X chromosome
d_act = 0.72; % degradation rate of free activator

K_n = 2.2; % concentration of bound SPEN at which activator synthesis rate is half max
n = 9.4;% Hill coefficient for SPEN supressing activator synthesis rate. 

a_x = 3.2/V; % xist synthesis rate from single X chromosome
d_x = 0.0076; % degradation rate of Xist 

m = 3.09; % Hill coefficient for bound SPEN reducing dissociation rate Xist
K_S = 2.06; % concentration of SPEN at which Xist dissociation is half max

K_a = 16.1; % concentration of activator at which Xist transcription is half max

k1 = 0.0084*V; % rate constant for Xist binding to DNA
k2 = 8.7; % maximum dissociation rate for Xist
k4 = 10.2; % dissoication rate for bound SPEN
k3 = 0.00026*V; % association rate for SPEN

sT = new_sT/V;  % total SPEN concentration
XbsT =  100/V; % concentration of Xist binding sites

N_S = round(sT/XbsT); % Number of SPEN that bind to one Xist

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