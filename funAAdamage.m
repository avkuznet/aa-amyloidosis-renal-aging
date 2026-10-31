function [Xi,BA,Q,Vfrac] = funAAdamage(p,B,t)
% Damage measures of Sections 2.2-2.3 evaluated from the steady-state
% oligomer concentration B, using analytical solutions (27), (41) and (46).
Q = B.*t;                                             % Eq. (46)
I = B*(p.T12I/p.theta12B).*(-expm1(-t*log(2)/p.T12I));% Eq. (27)
Vfrac = I*p.MWAA/(p.a21*p.rofibril);                  % Eq. (41)
Xi = Q + p.omega*Vfrac;                               % Eq. (47)
BA = p.tAAonset + t + p.c17*Xi;                       % Eq. (49)
end
