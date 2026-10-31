function out = funAAsolve(p,Stotfun,tf,nt)
% Integrates Eqs. (15), (19), (20) and (26) for a prescribed, possibly
% time-dependent total SAA concentration Stotfun(t) (uM) and returns the
% derived damage measures of Sections 2.2 and 2.3.
%
% p must contain: Kd, Htot, hS, Akidneys, Vkidneys, kcl, k1, k2, theta12B,
% T12S, T12A, T12B, T12I, MWAA, a21, rofibril, omega, c17, tAAonset.

Sfreefun = @(t) 0.5*(-(p.Kd+p.Htot-Stotfun(t)) + ...
    sqrt((p.Kd+p.Htot-Stotfun(t)).^2+4*p.Kd*Stotfun(t)));

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep',tf/2000);
sol = ode45(@(t,y) funCDCMCB_AA_tv(t,y,p.hS,p.Akidneys,p.Vkidneys,Sfreefun, ...
    p.kcl,p.T12S,p.k1,p.k2,p.T12A,p.T12B,p.theta12B,p.T12I),[0 tf],[0 0 0 0],options);

t = linspace(0,tf,nt);
y = deval(sol,t);

out.t = t;
out.S = y(1,:);
out.A = y(2,:);
out.B = y(3,:);
out.I = y(4,:);
out.Stot = Stotfun(t);
out.Sfree = Sfreefun(t);

out.Vfrac = out.I*p.MWAA/(p.a21*p.rofibril);   % Eq. (41)
out.Q     = cumtrapz(t,out.B);                 % Eq. (45)
out.Qtot  = out.Q + p.omega*out.Vfrac;         % Eq. (47), combined burden
out.biolage = p.tAAonset + t + p.c17*out.Qtot; % Eq. (49)
out.biolagemin = p.tAAonset + t + p.c17*out.Q; % Eq. (50), irreversible floor
end
