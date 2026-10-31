clc, clear, close all
tic

% Reproduces Fig. S17 of the paper. Renal biological age versus calendar age
% for three values of the half-deposition time theta_1/2,B.
%
% One run produces both panels of each figure: panel (a) with biologically
% relevant half-lives (functional protein degradation) and panel (b) with
% infinite half-lives (impaired protein degradation). Each figure window is
% named after the paper figure and panel it reproduces, e.g. 'Fig. S17a'.
% Model: Eqs. (15), (19), (20) and (26) with initial conditions (32);
% right-hand side in funCDCMCB_AA.m. Numerical (ode45) and analytical
% solutions are plotted together.

for panel = 'ab'
if panel == 'a'   % biologically relevant half-lives
  T12A=3.6e3; T12B=1.73e5; T12I=3.16e7; T12S=3.6e3;
  ttl='functional protein clearance';
else              % infinite half-lives
  T12A=1e14; T12B=1e14; T12I=1e14; T12S=1e14;
  ttl='impaired protein degradation pathways';
end
fprintf('\n=== panel (%s): %s ===\n',panel,ttl);

  theta12B=8.64e4;

  Vkidneys=4.04e-4;
  Akidneys=0.6;

  hS=5.55e-7;
  % k1=1e-6;
  k1=1e-4;
  % k2=1e-6;
  k2=1e-4;
  kcl=9.65e-5;
  % kcl=9.65e-3;

  MWAA=1.2e4;
  a21=1e21;
  rofibril=1.35e-12;

  % omegaS=0.05;
  % Htot=30;
  % Stot=8;
  %
  % Kd=omegaS*(Htot-Stot+omegaS*Stot)/(1-omegaS)

  Kd=1.1789;
  Stot=8;
  Htot=30;
  omega=1.4e10;
  c17=0.334;
  tAAonset=50*(365*24*3600);

  tf=1.89e8;

  %%%%%%%%%%%%%%%%
  theta12B=8.64e3;
  % theta12B=8.64e4;
  % theta12B=8.64e5;

  Sfree=0.5*(-(Kd+Htot-Stot)+sqrt((Kd+Htot-Stot).^2+4*Kd*Stot));
  omegaS=Sfree./Stot
  Hfree=Htot-Stot+Sfree;
  omegaH=Hfree./Htot

  options = odeset('RelTol',1e-10,'AbsTol',1e-10);
  tspan = [0 tf];
  y0 = [0 0 0 0];
  sol = ode45(@(t,CDMB) funCDCMCB_AA(t,CDMB,hS,Akidneys,Vkidneys,Sfree,kcl,T12S, ...
  k1,k2,T12A,T12B,theta12B,T12I), tspan, y0,...
      options);
  tnum = linspace(0,tf,250);
  CDMBnum = deval(sol,tnum);

  % y1=[S]
  % y2=[A]
  % y3=[B]
  % y4=[I]

  Skidney_1=CDMBnum(1,:);
  A_1=CDMBnum(2,:);
  B_1=CDMBnum(3,:);
  I_1=CDMBnum(4,:);

  Skidney_1anal=Akidneys*hS*Sfree*T12S/...
      (Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));

  Skidney=Skidney_1anal;

  B_1anal=k2.^(-1).*T12A.^(-1).*(T12B+theta12B).^(-1).*log(4).^(-1).*( ...
    k2.*kcl.*Skidney.*T12A.*T12B.*theta12B+(-1).*k1.*T12A.*( ...
    T12B+theta12B).*log(2)+(-1).*(T12B+theta12B).*log(2).^2+ ...
    T12B.*theta12B.*(k2.^2.*kcl.^2.*Skidney.^2.*T12A.^2+T12B.^( ...
    -2).*theta12B.^(-2).*(T12B+theta12B).^2.*log(2).^2.*(k1.* ...
    T12A+log(2)).^2+k2.*kcl.*Skidney.*T12A.*T12B.^(-1).* ...
    theta12B.^(-1).*(T12B+theta12B).*(k1.*T12A+(-1).*log(2)).* ...
    log(4)).^(1/2));

  A_1anal=T12A*(-B_1anal*(T12B+theta12B)/(T12B*theta12B)+kcl*Skidney/log(2));

  % I_1anal=tnum.*B_1anal*log(2)/theta12B;

  I_1anal=(1-exp(-tnum.*log(2)/T12I))*B_1anal*T12I/theta12B;

  VfibrilVkidney_1=I_1*MWAA/(a21*rofibril);
  VfibrilVkidney_1anal=I_1anal*MWAA/(a21*rofibril);

  Q_1(1)=0;
  for i29=2:250
      firsti29Componentst = tnum(1:i29);
      firsti29ComponentsBnum = B_1(1:i29);
  Q_1(i29)=trapz(firsti29Componentst,firsti29ComponentsBnum);
  end

  Q_1anal=B_1anal.*tnum;

  Qtot_1=Q_1+omega*VfibrilVkidney_1;
  Qtot_1anal=Q_1anal+omega*VfibrilVkidney_1anal;

  biolage_1=tAAonset+tnum+c17*Qtot_1;
  biolage_1anal=tAAonset+tnum+c17*Qtot_1anal;

  %%%%%%%%%%%%%%%%
  % theta12B=8.64e3;
  theta12B=8.64e4;
  % theta12B=8.64e5;

  Sfree=0.5*(-(Kd+Htot-Stot)+sqrt((Kd+Htot-Stot).^2+4*Kd*Stot));
  omegaS=Sfree./Stot
  Hfree=Htot-Stot+Sfree;
  omegaH=Hfree./Htot

  options = odeset('RelTol',1e-10,'AbsTol',1e-10);
  tspan = [0 tf];
  y0 = [0 0 0 0];
  sol = ode45(@(t,CDMB) funCDCMCB_AA(t,CDMB,hS,Akidneys,Vkidneys,Sfree,kcl,T12S, ...
  k1,k2,T12A,T12B,theta12B,T12I), tspan, y0,...
      options);
  tnum = linspace(0,tf,250);
  CDMBnum = deval(sol,tnum);

  Skidney_2=CDMBnum(1,:);
  A_2=CDMBnum(2,:);
  B_2=CDMBnum(3,:);
  I_2=CDMBnum(4,:);

  Skidney_2anal=Akidneys*hS*Sfree*T12S/...
      (Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));

  Skidney=Skidney_2anal;

  B_2anal=k2.^(-1).*T12A.^(-1).*(T12B+theta12B).^(-1).*log(4).^(-1).*( ...
    k2.*kcl.*Skidney.*T12A.*T12B.*theta12B+(-1).*k1.*T12A.*( ...
    T12B+theta12B).*log(2)+(-1).*(T12B+theta12B).*log(2).^2+ ...
    T12B.*theta12B.*(k2.^2.*kcl.^2.*Skidney.^2.*T12A.^2+T12B.^( ...
    -2).*theta12B.^(-2).*(T12B+theta12B).^2.*log(2).^2.*(k1.* ...
    T12A+log(2)).^2+k2.*kcl.*Skidney.*T12A.*T12B.^(-1).* ...
    theta12B.^(-1).*(T12B+theta12B).*(k1.*T12A+(-1).*log(2)).* ...
    log(4)).^(1/2));

  A_2anal=T12A*(-B_2anal*(T12B+theta12B)/(T12B*theta12B)+kcl*Skidney/log(2));

  % I_2anal=tnum.*B_2anal*log(2)/theta12B;

  I_2anal=(1-exp(-tnum.*log(2)/T12I))*B_2anal*T12I/theta12B;

  VfibrilVkidney_2=I_2*MWAA/(a21*rofibril);
  VfibrilVkidney_2anal=I_2anal*MWAA/(a21*rofibril);

  Q_2(1)=0;
  for i29=2:250
      firsti29Componentst = tnum(1:i29);
      firsti29ComponentsBnum = B_2(1:i29);
  Q_2(i29)=trapz(firsti29Componentst,firsti29ComponentsBnum);
  end

  Q_2anal=B_2anal.*tnum;

  Qtot_2=Q_2+omega*VfibrilVkidney_2;
  Qtot_2anal=Q_2anal+omega*VfibrilVkidney_2anal;

  biolage_2=tAAonset+tnum+c17*Qtot_2;
  biolage_2anal=tAAonset+tnum+c17*Qtot_2anal;

  %%%%%%%%%%%%%%%%
  % theta12B=8.64e3;
  % theta12B=8.64e4;
  theta12B=8.64e5;

  Sfree=0.5*(-(Kd+Htot-Stot)+sqrt((Kd+Htot-Stot).^2+4*Kd*Stot));
  omegaS=Sfree./Stot
  Hfree=Htot-Stot+Sfree;
  omegaH=Hfree./Htot

  options = odeset('RelTol',1e-10,'AbsTol',1e-10);
  tspan = [0 tf];
  y0 = [0 0 0 0];
  sol = ode45(@(t,CDMB) funCDCMCB_AA(t,CDMB,hS,Akidneys,Vkidneys,Sfree,kcl,T12S, ...
  k1,k2,T12A,T12B,theta12B,T12I), tspan, y0,...
      options);
  tnum = linspace(0,tf,250);
  CDMBnum = deval(sol,tnum);

  Skidney_3=CDMBnum(1,:);
  A_3=CDMBnum(2,:);
  B_3=CDMBnum(3,:);
  I_3=CDMBnum(4,:);

  Skidney_3anal=Akidneys*hS*Sfree*T12S/...
      (Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));

  Skidney=Skidney_3anal;

  B_3anal=k2.^(-1).*T12A.^(-1).*(T12B+theta12B).^(-1).*log(4).^(-1).*( ...
    k2.*kcl.*Skidney.*T12A.*T12B.*theta12B+(-1).*k1.*T12A.*( ...
    T12B+theta12B).*log(2)+(-1).*(T12B+theta12B).*log(2).^2+ ...
    T12B.*theta12B.*(k2.^2.*kcl.^2.*Skidney.^2.*T12A.^2+T12B.^( ...
    -2).*theta12B.^(-2).*(T12B+theta12B).^2.*log(2).^2.*(k1.* ...
    T12A+log(2)).^2+k2.*kcl.*Skidney.*T12A.*T12B.^(-1).* ...
    theta12B.^(-1).*(T12B+theta12B).*(k1.*T12A+(-1).*log(2)).* ...
    log(4)).^(1/2));

  A_3anal=T12A*(-B_3anal*(T12B+theta12B)/(T12B*theta12B)+kcl*Skidney/log(2));

  % I_3anal=tnum.*B_3anal*log(2)/theta12B;

  I_3anal=(1-exp(-tnum.*log(2)/T12I))*B_3anal*T12I/theta12B;

  VfibrilVkidney_3=I_3*MWAA/(a21*rofibril);
  VfibrilVkidney_3anal=I_3anal*MWAA/(a21*rofibril);

  Q_3(1)=0;
  for i29=2:250
      firsti29Componentst = tnum(1:i29);
      firsti29ComponentsBnum = B_3(1:i29);
  Q_3(i29)=trapz(firsti29Componentst,firsti29ComponentsBnum);
  end

  Q_3anal=B_3anal.*tnum;

  Qtot_3=Q_3+omega*VfibrilVkidney_3;
  Qtot_3anal=Q_3anal+omega*VfibrilVkidney_3anal;

  biolage_3=tAAonset+tnum+c17*Qtot_3;
  biolage_3anal=tAAonset+tnum+c17*Qtot_3anal;

  % y1=[S]kidney
  % y2=[A]
  % y3=[B]
  % y4=[I]

  figure('Name',['Fig. S17' panel],'NumberTitle','off')
  plot((tAAonset+tnum)/(365*24*3600),biolage_1/(365*24*3600),'-rx','LineWidth',1.5,'MarkerIndices',1:15:length(tnum));
  hold on
  plot((tAAonset+tnum)/(365*24*3600),biolage_1anal/(365*24*3600),'-c*','LineWidth',1.5,'MarkerIndices',1:17:length(tnum));
  plot((tAAonset+tnum)/(365*24*3600),biolage_2/(365*24*3600),'--bsquare','LineWidth',1.5,'MarkerIndices',1:20:length(tnum));
  plot((tAAonset+tnum)/(365*24*3600),biolage_2anal/(365*24*3600),'--m>','LineWidth',1.5,'MarkerIndices',1:22:length(tnum));
  plot((tAAonset+tnum)/(365*24*3600),biolage_3/(365*24*3600),'-.g^','LineWidth',1.5,'MarkerIndices',1:25:length(tnum));
  plot((tAAonset+tnum)/(365*24*3600),biolage_3anal/(365*24*3600),'-.ko','LineWidth',1.5,'MarkerIndices',1:27:length(tnum));
  plot((tAAonset+tnum)/(365*24*3600),100.*tnum./tnum,':r^','LineWidth',1.5,'MarkerIndices',1:29:length(tnum));
  hold off
  title(ttl)
  xlabel('calendar age (years)');
  ylabel('biological age (years)');
  legend('\theta_{1/2,B}=8.64\times10^3 s numerical','\theta_{1/2,B}=8.64\times10^3 s analytical',...
      '\theta_{1/2,B}=8.64\times10^4 s numerical','\theta_{1/2,B}=8.64\times10^4 s analytical',...
      '\theta_{1/2,B}=8.64\times10^5 s numerical','\theta_{1/2,B}=8.64\times10^5 s analytical',...
      '100 years','location','best')
  % ylim([50 100])
  ylim([50 inf])

  Q_1(end)
  VfibrilVkidney_1(end)
  Qtot_1(end)
  biolage_1(end)/(365*24*3600)
end

toc
