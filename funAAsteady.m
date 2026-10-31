function s = funAAsteady(p,Stotv)
% Closed-form steady state of Eqs. (15),(19),(20) for a constant total SAA
% concentration, written in terms of the two dimensionless groups
%   kappa1 = k1*T12A/ln2,  kappa2 = k2*Bstar*T12A/ln2
% with Bstar the supply-limited oligomer concentration.

s.Sfree = 0.5*(-(p.Kd+p.Htot-Stotv)+sqrt((p.Kd+p.Htot-Stotv).^2+4*p.Kd*Stotv));
s.S = p.Akidneys*p.hS*s.Sfree*p.T12S./ ...
      (p.Akidneys*p.hS*p.T12S+p.kcl*p.T12S*p.Vkidneys+p.Vkidneys*log(2));
s.Bstar = p.kcl*s.S*p.T12B*p.theta12B/(log(2)*(p.T12B+p.theta12B));
s.kappa1 = p.k1*p.T12A/log(2);
s.kappa2 = p.k2*s.Bstar*p.T12A/log(2);
s.b = (-(1+s.kappa1-s.kappa2)+sqrt((1+s.kappa1-s.kappa2).^2+4*s.kappa1.*s.kappa2)) ...
      ./(2*s.kappa2);
s.B = s.b.*s.Bstar;
s.A = (p.T12A/log(2))*(p.kcl*s.S-s.B*log(2)*(p.T12B+p.theta12B)/(p.T12B*p.theta12B));
end
