function p = xci_scale_volume(p, V)
%XCI_SCALE_VOLUME  Rescale the standard parameter vector to concentration units.
%
%   p = XCI_SCALE_VOLUME(p, V) converts the quantity-based parameters to
%   concentration units for a nuclear volume V (used by Fig3B and Fig3C).
%   Intensive parameters are unchanged; extensive ones scale with V:
%       a_act -> a_act/V,  a_x -> a_x/V,  k1 -> k1*V,  k3 -> k3*V,  XbsT -> XbsT/V
%
%   sT is NOT scaled here because different scripts set it differently
%   (Fig3B uses sT/V, Fig3C uses new_sT/V); the caller overwrites p(15).

p(1)  = p(1)/V;    % a_act
p(5)  = p(5)/V;    % a_x
p(10) = p(10)*V;   % k1
p(13) = p(13)*V;   % k3
p(14) = p(14)/V;   % XbsT
end
