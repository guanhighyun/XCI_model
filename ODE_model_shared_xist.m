function dy = ODE_model_shared_xist(t, y, p, N_c)
%ODE_MODEL_SHARED_XIST  XCI model variant with a single shared free-Xist pool.
%
%   Same biology as ODE_MODEL, but all chromosomes draw from ONE pooled free
%   Xist species rather than a per-chromosome free-Xist pool (used by FigS6
%   and FigS7). One activator and one free-Xist pool are shared; each
%   chromosome keeps its own bound Xist and bound SPEN.
%
%   State vector y (length 2 + 2*N_c):
%       y(1)      = act        free activator (shared)
%       y(2)      = xf         free Xist (shared)
%     for chromosome i = 1..N_c:
%       y(2*i+1)  = x_i^b      bound Xist
%       y(2*i+2)  = s_i^b      bound SPEN
%
%   p is the 15-element parameter vector from XCI_PARAMS. At N_c = 1 this is
%   identical to the standard ODE_MODEL (a single chromosome shares its own
%   free Xist).

a_act = p(1);  d_act = p(2);  K_n = p(3);  n = p(4);
a_x   = p(5);  d_x   = p(6);  m   = p(7);  K_S = p(8);  K_a = p(9);
k1    = p(10); k2    = p(11); k4  = p(12); k3  = p(13);
XbsT  = p(14); sT    = p(15);

N_S = round(sT/XbsT);

act = y(1);
xf  = y(2);

% Total bound SPEN across chromosomes -> free SPEN by global conservation.
Sb_total = 0;
for i = 1:N_c
    Sb_total = Sb_total + y(2*i+2);
end
free_SPEN = sT - Sb_total;

dy = zeros(2 + 2*N_c, 1);

% Shared activator: basal synthesis repressed by each chromosome's bound SPEN.
dact = -d_act*act;
for i = 1:N_c
    sib = y(2*i+2);
    dact = dact + a_act/(1 + (sib/K_n)^n);
end
dy(1) = dact;

% Shared free Xist: synthesised/degraded once, exchanged with every chromosome.
dxf = a_x*act/(K_a + act) - d_x*xf;
for i = 1:N_c
    xib = y(2*i+1);
    sib = y(2*i+2);
    dxf = dxf - k1*(XbsT - xib)*xf + k2*xib/(1 + (sib/K_S)^m);
end
dy(2) = dxf;

% Per-chromosome bound Xist and bound SPEN.
for i = 1:N_c
    xib = y(2*i+1);
    sib = y(2*i+2);
    dy(2*i+1) = k1*(XbsT - xib)*xf - k2*xib/(1 + (sib/K_S)^m);   % bound Xist
    dy(2*i+2) = k3*free_SPEN*(N_S*xib - sib) - k4*sib;           % bound SPEN
end
end
