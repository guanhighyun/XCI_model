function dy = ODE_model(t, y, p, N_c)
%ODE_MODEL  Right-hand side of the XCI model for N_c X chromosomes.
%
%   dy = ODE_MODEL(t, y, p, N_c) returns the time derivatives of the state
%   vector y for a cell with N_c X chromosomes. One activator pool is shared
%   across all chromosomes; each chromosome has its own free Xist, bound
%   Xist and bound SPEN.
%
%   State vector y (length 1 + 3*N_c):
%       y(1)      = act        free activator (shared)
%     for chromosome i = 1..N_c:
%       y(3*i-1)  = x_i^f      free Xist
%       y(3*i)    = x_i^b      bound Xist
%       y(3*i+1)  = s_i^b      bound SPEN
%
%   p is the 15-element parameter vector from XCI_PARAMS (see that function
%   for the element order). N_c = 1,2,3,4 reproduce the XY/XX/XXX/tetraploid
%   cases exactly.
%
%   Setting N_c = 1,2,3 reproduces the Python DE1/DE2/DE3 right-hand sides.

a_act = p(1);  d_act = p(2);  K_n = p(3);  n = p(4);
a_x   = p(5);  d_x   = p(6);  m   = p(7);  K_S = p(8);  K_a = p(9);
k1    = p(10); k2    = p(11); k4  = p(12); k3  = p(13);
XbsT  = p(14); sT    = p(15);

% Number of SPEN that can bind one Xist (each chromosome can recruit all SPEN).
N_S = round(sT/XbsT);

act = y(1);

% Total bound SPEN across chromosomes -> free SPEN by global conservation.
Sb_total = 0;
for i = 1:N_c
    Sb_total = Sb_total + y(3*i+1);
end
free_SPEN = sT - Sb_total;

dy = zeros(1 + 3*N_c, 1);

% Shared activator: basal synthesis repressed by each chromosome's bound SPEN.
dact = -d_act*act;
for i = 1:N_c
    sib = y(3*i+1);
    dact = dact + a_act/(1 + (sib/K_n)^n);
end
dy(1) = dact;

% Per-chromosome Xist and SPEN dynamics.
for i = 1:N_c
    xif = y(3*i-1);
    xib = y(3*i);
    sib = y(3*i+1);
    bind   = k1*(XbsT - xib)*xif;              % Xist binding to DNA
    unbind = k2*xib/(1 + (sib/K_S)^m);          % dissociation, slowed by bound SPEN
    dy(3*i-1) = a_x*act/(K_a + act) - d_x*xif - bind + unbind;   % free Xist
    dy(3*i)   = bind - unbind;                                   % bound Xist
    dy(3*i+1) = k3*free_SPEN*(N_S*xib - sib) - k4*sib;           % bound SPEN
end
end
