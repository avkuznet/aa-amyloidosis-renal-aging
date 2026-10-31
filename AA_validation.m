clc, clear, close all
tic
% Reproduces Fig. 11 of the paper.
% Comparison of model-predicted progression rates with the clinical
% SAA-outcome relationships reported by Gillmore et al. (2001) and
% Lachmann et al. (2007). Uses analytical solutions (16), (25), (27), (46).

AA_params

p.Kd=Kd; p.Htot=Htot; p.hS=hS; p.Akidneys=Akidneys; p.Vkidneys=Vkidneys;
p.kcl=kcl; p.k1=k1; p.k2=k2; p.theta12B=theta12B;
p.MWAA=MWAA; p.a21=a21; p.rofibril=rofibril;
p.omega=omega; p.c17=c17; p.tAAonset=tAAonset;

mgL2uM=1e3/MWAA;                 % mg/L -> uM for a 12 kDa protein
SmgL=logspace(log10(2),log10(200),200);
SuM=SmgL*mgL2uM;

%%%%%%%%%%%%%%%% verification of the analytical route against ode45
p.T12S=T12S_f; p.T12A=T12A_f; p.T12B=T12B_f; p.T12I=T12I_f;
for Sv=[2 8]
  ss=funAAsteady(p,Sv);
  [Xi_a,BA_a]=funAAdamage(p,ss.B,6*yr);
  r=funAAsolve(p,@(t) Sv*ones(size(t)),6*yr,3000);
  fprintf(['check [S]tot=%4.1f uM : Xi analytical=%11.5e numerical=%11.5e ' ...
      '(rel.diff %7.2e), biol.age %8.4f vs %8.4f yr\n'], ...
      Sv,Xi_a,r.Qtot(end),abs(Xi_a-r.Qtot(end))/r.Qtot(end),BA_a/yr,r.biolage(end)/yr);
end

%%%%%%%%%%%%%%%% time to the calibrated terminal-damage threshold, BA = 100 yr
t100=nan(2,numel(SuM));
for scen=1:2
  if scen==1
    p.T12S=T12S_f; p.T12A=T12A_f; p.T12B=T12B_f; p.T12I=T12I_f;
  else
    p.T12S=T12S_i; p.T12A=T12A_i; p.T12B=T12B_i; p.T12I=T12I_i;
  end
  for i=1:numel(SuM)
    ss=funAAsteady(p,SuM(i));
    f=@(t) p.tAAonset+t+p.c17*(ss.B*t+p.omega*ss.B*(p.T12I/p.theta12B)* ...
        (-expm1(-t*log(2)/p.T12I))*p.MWAA/(p.a21*p.rofibril))-100*yr;
    if f(200*yr)<0, continue, end
    t100(scen,i)=fzero(f,[1e-6 200*yr])/yr;
  end
end

fprintf('\n SAA(mg/L)  [S]tot(uM)   t100 functional(yr)   t100 impaired(yr)\n');
for v=[2 5 10 20 30 50 100 200]
  i=find(SmgL>=v,1);
  fprintf('%9.0f %11.3f %20.2f %19.2f\n',SmgL(i),SuM(i),t100(1,i),t100(2,i));
end

rate=1./t100;
ref=10; iref=[interp1(SmgL,rate(1,:),ref) interp1(SmgL,rate(2,:),ref)];
fprintf('\nfold increase in the rate of renal biological ageing per doubling of SAA\n');
nm={'functional clearance','impaired degradation'};
for scen=1:2
  lg=log2(SmgL); lr=log(rate(scen,:)); ok=isfinite(lr);
  for band=[10 20;20 40;40 80;10 100]'
    m=ok & SmgL>=band(1) & SmgL<=band(2);
    pf=polyfit(lg(m),lr(m),1);
    fprintf('  %-22s %4.0f-%4.0f mg/L : %6.3f per doubling\n',nm{scen},band(1),band(2),exp(pf(1)));
  end
end
fprintf('  reported by Lachmann et al. (2007)            : 1.24 per doubling (ESRD risk)\n');

figure(1)
loglog(SmgL,t100(1,:),'-r','LineWidth',1.5);
hold on
loglog(SmgL,t100(2,:),'--b','LineWidth',1.5);
yl=[min(t100(:)) max(t100(:))];
loglog([10 10],yl,':k','LineWidth',1.5);
loglog([50 50],yl,'-.m','LineWidth',1.5);
hold off
xlabel('median circulating SAA (mg/L)');
ylabel('time to BA = 100 years (years)');
legend('functional protein clearance','impaired protein degradation pathways',...
    'SAA = 10 mg/L','SAA = 50 mg/L','location','best');

figure(2)
semilogx(SmgL,rate(1,:)/iref(1),'-r','LineWidth',1.5);
hold on
semilogx(SmgL,rate(2,:)/iref(2),'--b','LineWidth',1.5);
semilogx(SmgL,1.24.^log2(SmgL/ref),':k','LineWidth',2);
hold off
xlabel('median circulating SAA (mg/L)');
ylabel('rate of renal biological ageing relative to SAA = 10 mg/L');
legend('functional protein clearance','impaired protein degradation pathways',...
    'Lachmann et al. (2007), 1.24 per doubling','location','best');

set(figure(1),'Name','Fig. 11a','NumberTitle','off');
savefig_emf(figure(1),'fig11a');
set(figure(2),'Name','Fig. 11b','NumberTitle','off');
savefig_emf(figure(2),'fig11b');

toc
