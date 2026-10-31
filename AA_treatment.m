clc, clear, close all
tic
% Reproduces Fig. 10 of the paper.
% Treatment-response simulation: suppression of circulating SAA to the
% target of ref. (Gillmore et al., 2001) (<10 mg/L = 0.83 uM) started at
% different times after disease onset.

AA_params

p.Kd=Kd; p.Htot=Htot; p.hS=hS; p.Akidneys=Akidneys; p.Vkidneys=Vkidneys;
p.kcl=kcl; p.k1=k1; p.k2=k2; p.theta12B=theta12B;
p.MWAA=MWAA; p.a21=a21; p.rofibril=rofibril;
p.omega=omega; p.c17=c17; p.tAAonset=tAAonset;

tf=10*yr; nt=4000;
Shigh=8; Starget=0.83;       % uM; 0.83 uM = 10 mg/L
w=30*24*3600;                % 1 month to reach the SAA target on therapy
trx=[1 3 5]*yr;
sty={'--bdiamond','-.g^',':mv'};

for scen=1:2
  if scen==1
    p.T12S=T12S_f; p.T12A=T12A_f; p.T12B=T12B_f; p.T12I=T12I_f;
    ttl='functional protein clearance';
  else
    p.T12S=T12S_i; p.T12A=T12A_i; p.T12B=T12B_i; p.T12I=T12I_i;
    ttl='impaired protein degradation pathways';
  end
  fprintf('\n--- %s ---\n',ttl);

  figure(scen); hold on
  r0=funAAsolve(p,@(t) Shigh*ones(size(t)),tf,nt);
  plot(r0.t/yr,r0.biolage/yr,'-rs','LineWidth',1.5,'MarkerIndices',1:300:nt);
  fprintf('no therapy                 biol.age(10yr)=%8.2f yr\n',r0.biolage(end)/yr);
  lab={'no therapy'};
  for j=1:numel(trx)
    t0=trx(j);
    Sfun=@(t) Starget+(Shigh-Starget)./(1+exp((t-t0)/w));
    r=funAAsolve(p,Sfun,tf,nt);
    plot(r.t/yr,r.biolage/yr,sty{j},'LineWidth',1.5,'MarkerIndices',1:300+40*j:nt);
    lab{end+1}=sprintf('SAA suppressed from t = %g yr',t0/yr);
    fprintf(['therapy at %g yr            biol.age(10yr)=%8.2f yr   ' ...
        'biol.age at start of therapy=%8.2f yr\n'], ...
        t0/yr,r.biolage(end)/yr,interp1(r.t,r.biolage,t0)/yr);
  end
  plot(r0.t/yr,100*ones(size(r0.t)),':k','LineWidth',1.5);
  lab{end+1}='calibrated threshold (BA = 100 years)';
  hold off
  xlabel('t (years)'); ylabel('biological age (years)');
  title(ttl); legend(lab,'location','best');
end

set(figure(1),'Name','Fig. 10a','NumberTitle','off');
savefig_emf(figure(1),'fig10a');
set(figure(2),'Name','Fig. 10b','NumberTitle','off');
savefig_emf(figure(2),'fig10b');

toc
