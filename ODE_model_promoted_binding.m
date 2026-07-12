function dy = ODE_model_promoted_binding(t, y, p, N_c, suppress_dissoc)
%ODE_MODEL_PROMOTED_BINDING  XCI model variant where bound SPEN promotes Xist binding.
%
%   Same state layout as the standard ODE_MODEL (per-chromosome free Xist),
%   but the Xist binding rate is increased by bound SPEN via an extra k5 term
%   (used by FigS3 and FigS4).
%
%   State vector y (length 1 + 3*N_c):
%       y(1)      = act        free activator (shared)
%     for chromosome i = 1..N_c:
%       y(3*i-1)  = x_i^f      free Xist
%       y(3*i)    = x_i^b      bound Xist
%       y(3*i+1)  = s_i^b      bound SPEN
%
%   p is the 16-element parameter vector from XCI_PARAMS('figS3' | 'figS4');
%   element 16 is k5 (SPEN-promoted binding rate).
%
%   suppress_dissoc selects the Xist dissociation form:
%       false -> k2*x_i^b                          (FigS3)
%       true  -> k2*x_i^b/(1 + (s_i^b/K_S)^m)       (FigS4)

a_act = p(1);  d_act = p(2);  K_n = p(3);  n = p(4);
a_x   = p(5);  d_x   = p(6);  m   = p(7);  K_S = p(8);  K_a = p(9);
k1    = p(10); k2    = p(11); k4  = p(12); k3  = p(13);
XbsT  = p(14); sT    = p(15); k5  = p(16);

N_S = round(sT/XbsT);

act = y(1);

Sb_total = 0;
for i = 1:N_c
    Sb_total = Sb_total + y(3*i+1);
end
free_SPEN = sT - Sb_total;

dy = zeros(1 + 3*N_c, 1);

dact = -d_act*act;
for i = 1:N_c
    sib = y(3*i+1);
    dact = dact + a_act/(1 + (sib/K_n)^n);
end
dy(1) = dact;

for i = 1:N_c
    xif = y(3*i-1);
    xib = y(3*i);
    sib = y(3*i+1);
    bind = (XbsT - xib)*xif*(k1 + k5*(sib^m)/(sib^m + K_S^m));   % SPEN-promoted binding
    if suppress_dissoc
        dissoc = k2*xib/(1 + (sib/K_S)^m);
    else
        dissoc = k2*xib;
    end
    dy(3*i-1) = a_x*act/(K_a + act) - d_x*xif - bind + dissoc;   % free Xist
    dy(3*i)   = bind - dissoc;                                   % bound Xist
    dy(3*i+1) = k3*free_SPEN*(N_S*xib - sib) - k4*sib;           % bound SPEN
end
end
