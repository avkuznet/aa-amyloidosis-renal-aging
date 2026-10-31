clc, clear, close all
tic
% Reproduces Figs. 9 and S27 of the paper.
% Intermittent inflammatory flares versus sustained inflammation at the same
% time-averaged SAA concentration (Fig. 9), and a single inflammatory episode
% placed early or late in the disease course (Fig. S27). Tests whether the
% path dependence of Xi_comb is a dynamical consequence of the model or
% merely a restatement of its definition. Panels (a) of both figures show
% the prescribed total SAA concentration; Figs. 9b-d show the free
% (HDL-unbound) SAA concentration and the renal oligomer concentration [B]
% that result, and Figs. 9e,f and S27b,c the combined renal damage.

AA_params

p.Kd=Kd; p.Htot=Htot; p.hS=hS; p.Akidneys=Akidneys; p.Vkidneys=Vkidneys;
p.kcl=kcl; p.k1=k1; p.k2=k2; p.theta12B=theta12B;
p.MWAA=MWAA; p.a21=a21; p.rofibril=rofibril;
p.omega=omega; p.c17=c17; p.tAAonset=tAAonset;

tf=6*yr; nt=4000;
Smean=4;                 % time-averaged total SAA, uM
period=yr;               % one flare per year
Slo=0.5;                 % inter-flare (controlled) SAA, uM
w=2*24*3600;             % 2-day rise/fall time of the acute-phase response
ton=0.15*period;

duty=[1 0.25 0.10];
lab={'sustained ([S]_{tot}=4 \muM)','flares, 25% duty cycle','flares, 10% duty cycle'};
sty={'-rs','--bdiamond','-.g^'};

for scen=1:2
  if scen==1
    p.T12S=T12S_f; p.T12A=T12A_f; p.T12B=T12B_f; p.T12I=T12I_f;
    ttl='functional protein clearance';
  else
    p.T12S=T12S_i; p.T12A=T12A_i; p.T12B=T12B_i; p.T12I=T12I_i;
    ttl='impaired protein degradation pathways';
  end

  figure(scen); hold on
  fprintf('\n--- %s ---\n',ttl);
  base=NaN;
  for j=1:numel(duty)
    d=duty(j);
    if d==1
      Sfun=@(t) Smean*ones(size(t));
    else
      Shi0=(Smean-(1-d)*Slo)/d;
      tt=linspace(0,period,20001);
      m0=trapz(tt,funflare(tt,Slo,Shi0,period,d,ton,w))/period;
      Shi=Slo+(Shi0-Slo)*(Smean-Slo)/(m0-Slo);   % rescale to hit the target mean
      Sfun=@(t) funflare(t,Slo,Shi,period,d,ton,w);
    end
    Sin{j}=Sfun;
    r=funAAsolve(p,Sfun,tf,nt);
    R{scen,j}=r;
    plot(r.t/yr,r.Qtot,sty{j},'LineWidth',1.5,'MarkerIndices',1:300+40*j:nt);
    if d==1, Shi=Smean; end
    fprintf(['duty=%5.2f  peak[S]tot=%7.3f uM  realized mean=%7.4f uM  ' ...
        'Xi_comb(6yr)=%10.4e uM s  biol.age(6yr)=%7.3f yr\n'], ...
        d,Shi,trapz(r.t,r.Stot)/tf,r.Qtot(end),r.biolage(end)/yr);
    if j==1, base=r.Qtot(end); else
      fprintf('            ratio to sustained = %6.3f\n',r.Qtot(end)/base); end
  end
  hold off
  xlabel('t (years)'); ylabel('\Xi_{comb} (\muM s)');
  title(ttl); legend(lab,'location','best');
end

% prescribed total SAA of the three histories (same for both clearance cases)
tp=linspace(0,tf,6001);
figure(5); hold on
for j=1:numel(duty)
  plot(tp/yr,Sin{j}(tp),sty{j},'LineWidth',1.5,'MarkerIndices',1:300+40*j:numel(tp));
end
hold off
xlabel('t (years)'); ylabel('[S]_{tot} (\muM)');
ylim([0 60]);
title('prescribed total SAA'); legend(lab,'location','north');

% what the kidney receives: free (HDL-unbound) SAA, Eq. (9), and the
% resulting renal oligomer concentration [B]
fprintf('\n--- time averages over six years ---\n');
figure(7); hold on
for j=1:numel(duty)
  r=R{1,j};
  semilogy(r.t/yr,r.Sfree,sty{j},'LineWidth',1.5,'MarkerIndices',1:300+40*j:nt);
  fprintf('%-30s mean[S]tot=%7.4f  mean[S]free=%7.4f uM\n', ...
      lab{j},trapz(r.t,r.Stot)/tf,trapz(r.t,r.Sfree)/tf);
end
hold off; set(gca,'YScale','log'); ylim([1e-2 1e3]);
xlabel('t (years)'); ylabel('[S]_{free} (\muM)');
title('free SAA in plasma'); legend(lab,'location','north');
for scen=1:2
  figure(7+scen); hold on
  for j=1:numel(duty)
    r=R{scen,j};
    semilogy(r.t/yr,r.B,sty{j},'LineWidth',1.5,'MarkerIndices',1:300+40*j:nt);
    fprintf('scenario %d  %-30s mean[B]=%10.4e uM\n',scen,lab{j},trapz(r.t,r.B)/tf);
  end
  hold off; set(gca,'YScale','log'); ylim([1e-2 1e4]);
  xlabel('t (years)'); ylabel('[B] (\muM)');
  if scen==1, title('functional protein clearance');
  else, title('impaired protein degradation pathways'); end
  legend(lab,'location','north');
end

%%%%%%%%%%%%%%%% early versus late inflammatory episode of identical magnitude
fprintf('\n--- early vs late inflammatory episode (identical dose) ---\n');
Shi=16; dur=yr; onsets=[0.5*yr 4.5*yr];
onlab={'early episode (0.5-1.5 yr)','late episode (4.5-5.5 yr)'};
onsty={'-rs','--bdiamond'};
for scen=1:2
  if scen==1
    p.T12S=T12S_f; p.T12A=T12A_f; p.T12B=T12B_f; p.T12I=T12I_f;
    ttl='functional protein clearance';
  else
    p.T12S=T12S_i; p.T12A=T12A_i; p.T12B=T12B_i; p.T12I=T12I_i;
    ttl='impaired protein degradation pathways';
  end
  figure(2+scen); hold on
  fprintf('  [%s]\n',ttl);
  for j=1:2
    t0=onsets(j);
    Sfun=@(t) funflare(t-t0+0.15*1e9,Slo,Shi,1e9,dur/1e9,0.15*1e9,w);
    Sep{j}=Sfun;
    r=funAAsolve(p,Sfun,tf,nt);
    plot(r.t/yr,r.Qtot,onsty{j},'LineWidth',1.5,'MarkerIndices',1:300+60*j:nt);
    fprintf(['  %-26s mean[S]tot=%7.4f  Xi_comb(6yr)=%10.4e  ' ...
        'biol.age(6yr)=%7.3f yr  V_fibril/V_kidney(6yr)=%7.4f\n'], ...
        onlab{j},trapz(r.t,r.Stot)/tf,r.Qtot(end),r.biolage(end)/yr,r.Vfrac(end));
  end
  hold off
  xlabel('t (years)'); ylabel('\Xi_{comb} (\muM s)');
  title(ttl); legend(onlab,'location','best');
end

figure(6); hold on
for j=1:2
  plot(tp/yr,Sep{j}(tp),onsty{j},'LineWidth',1.5,'MarkerIndices',1:300+60*j:numel(tp));
end
hold off
xlabel('t (years)'); ylabel('[S]_{tot} (\muM)');
title('prescribed total SAA'); legend(onlab,'location','best');

% panel order follows the causal chain: input, free SAA, oligomers, damage
set(figure(5),'Name','Fig. 9a','NumberTitle','off');
savefig_emf(figure(5),'fig9a');
set(figure(7),'Name','Fig. 9b','NumberTitle','off');
savefig_emf(figure(7),'fig9b');
set(figure(8),'Name','Fig. 9c','NumberTitle','off');
savefig_emf(figure(8),'fig9c');
set(figure(9),'Name','Fig. 9d','NumberTitle','off');
savefig_emf(figure(9),'fig9d');
set(figure(1),'Name','Fig. 9e','NumberTitle','off');
savefig_emf(figure(1),'fig9e');
set(figure(2),'Name','Fig. 9f','NumberTitle','off');
savefig_emf(figure(2),'fig9f');
set(figure(6),'Name','Fig. S27a','NumberTitle','off');
savefig_emf(figure(6),'figS27a');
set(figure(3),'Name','Fig. S27b','NumberTitle','off');
savefig_emf(figure(3),'figS27b');
set(figure(4),'Name','Fig. S27c','NumberTitle','off');
savefig_emf(figure(4),'figS27c');

toc
