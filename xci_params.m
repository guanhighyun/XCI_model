function p = xci_params(name)
%XCI_PARAMS  Canonical model parameter vector for a named fitted set.
%
%   p = XCI_PARAMS(NAME) returns the parameter vector used by the XCI model
%   functions (ODE_model, ODE_model_shared_xist, ODE_model_promoted_binding).
%
%   Element order (15 elements), matching the Python EA convention:
%       1  a_act   activator synthesis rate from a single X chromosome
%       2  d_act   degradation rate of free activator
%       3  K_n     bound SPEN at which activator synthesis is half-max
%       4  n       Hill coefficient, SPEN suppressing activator synthesis
%       5  a_x     Xist synthesis rate from a single X chromosome
%       6  d_x     degradation rate of Xist
%       7  m       Hill coefficient, bound SPEN reducing Xist dissociation
%       8  K_S     SPEN at which Xist dissociation is half-max
%       9  K_a     activator at which Xist transcription is half-max
%      10  k1      rate constant for Xist binding to DNA
%      11  k2      maximum dissociation rate for Xist
%      12  k4      dissociation rate for bound SPEN
%      13  k3      association rate for SPEN     (NOTE: k4 precedes k3)
%      14  XbsT    quantity of Xist binding sites (per chromosome)
%      15  sT      total SPEN quantity (global; may be overwritten by caller)
%
%   The 'figS3' and 'figS4' sets append a 16th element:
%      16  k5      rate at which bound SPEN promotes Xist binding
%
%   sT is returned at the default used by the corresponding figure; scripts
%   that sweep or rescale sT overwrite p(15) after calling this function.

switch name
    case 'activator_inhibition'   % Fig1, Fig3 (base set), FigS1
        p = [1.1 0.72 2.2 9.4 3.2 0.0076 3.09 2.06 16.1 0.0084 8.7 10.2 0.00026 100 1000];
    case 'fig2'                   % Fig2, FigS2
        p = [12.8 0.13 7.13 1.0 9.0 0.031 2.0 4.3 374.6 0.0020 4.54 0.11 8.26 100 70];
    case 'spen_depletion'         % Fig4
        p = [66.8 0.048 1.0 10.0 5.9 0.0041 3.85 7.40 76.5 0.0024 4.21 8.78 0.014 100 250];
    case 'figS6'                  % FigS6 (shared free-Xist variant)
        p = [1.2 4.3 1.0 10.0 7.4 0.014 2.58 1.23 2.19 0.00084 0.27 0.36 0.0015 100 52];
    case 'figS7'                  % FigS7 (shared free-Xist variant)
        % NOTE: the source scripts label k4/k3 with swapped comments; the
        % numeric values below match the code's actual use (k3 = association).
        p = [0.707 0.028 57.8 4.39 3.9 0.027 2.58 26.8 351 0.0018 1.42 0.10 9.9 101 676];
    case 'figS3'                  % k5 promoted-binding, plain dissociation (16 elems)
        p = [1.89 3.43 2.85 7.8 10.0 0.029 3.6 5.78 4.15 0.011 8.57 8.39 0.00019 101 676 1.28];
    case 'figS4'                  % k5 promoted-binding, SPEN-suppressed dissociation (16 elems)
        p = [37.8 0.12 2.23 9.8 6.01 0.028 9.6 19.9 31.1 0.00010 9.28 0.11 5.91 124 1303 2.57];
    otherwise
        error('xci_params:unknownSet', 'Unknown parameter set: %s', name);
end
end
