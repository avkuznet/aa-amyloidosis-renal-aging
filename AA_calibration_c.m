clc, clear, close all
tic
% Reproduces Fig. S28 of the paper.
% Robustness of the renal biological-age scale to the calibration of c.
% c is fixed by back-solving Eq. (49) from an assumed terminal scenario, so
% this script quantifies how much of the model output that choice controls.

AA_params

p.Kd=Kd; p.Htot=Htot; p.hS=hS; p.Akidneys=Akidneys; p.Vkidneys=Vkidneys;
p.kcl=kcl; p.k1=k1; p.k2=k2; p.theta12B=theta12B;
p.MWAA=MWAA; p.a21=a21; p.rofibril=rofibril;
p.omega=omega; p.c17=c17; p.tAAonset=tAAonset;
p.T12S=T12S_i; p.T12A=T12A_i; p.T12B=T12B_i; p.T12I=T12I_i;

nt=3000;

% ---- consistency check: the published calibration must reproduce 100 yr ----
r=funAAsolve(p,@(t) Stot*ones(size(t)),6*yr,nt);
fprintf('check: biological age at 6 yr with published c = %.3f yr (target 100)\n',...
    r.biolage(end)/yr);
fprintf('       Xi_comb(6 yr) = %.4e uM s,  V_fibril/V_kidney = %.4f\n\n',...
    r.Qtot(end),r.Vfrac(end));

% ---- c implied by the reported 4-8 yr survival range ----
tterm=(4:8)*yr;
cv=zeros(size(tterm));
fprintf('assumed terminal survival   implied c (uM^-1)   ratio to published c\n');
for i=1:numel(tterm)
  ri=funAAsolve(p,@(t) Stot*ones(size(t)),tterm(i),nt);
  cv(i)=(100*yr-tAAonset-tterm(i))/ri.Qtot(end);
  fprintf('%17.0f yr %18.4f %21.3f\n',tterm(i)/yr,cv(i),cv(i)/c17);
end

% ---- effect of c on the dimensionless sensitivity coefficients ----
% evaluated at the Table 3 baseline ([H]tot = 30 uM); note that c itself is
% calibrated at [H]tot = 15 uM, as described in Section 2.3.
p.Htot=30;
delta=0.001;
pars={'Htot','kcl','hS','theta12B','Kd','Stot','k1','k2'};
tev=6*yr;
Sens=zeros(numel(cv),numel(pars));

for i=1:numel(cv)
  p.c17=cv(i);
  base=funAAsolve(p,@(t) Stot*ones(size(t)),tev,nt);
  BA0=base.biolage(end);
  for j=1:numel(pars)
    q=p; Sl=Stot;
    switch pars{j}
      case 'Htot',     q.Htot=p.Htot*(1+delta);
      case 'kcl',      q.kcl=p.kcl*(1+delta);
      case 'hS',       q.hS=p.hS*(1+delta);
      case 'theta12B', q.theta12B=p.theta12B*(1+delta);
      case 'Kd',       q.Kd=p.Kd*(1+delta);
      case 'Stot',     Sl=Stot*(1+delta);
      case 'k1',       q.k1=p.k1*(1+delta);
      case 'k2',       q.k2=p.k2*(1+delta);
    end
    rp=funAAsolve(q,@(t) Sl*ones(size(t)),tev,nt);
    Sens(i,j)=(1/BA0)*((rp.biolage(end)-BA0)/delta);
  end
end
p.c17=c17;

fprintf('\ndimensionless sensitivity of biological age at 6 yr\n');
fprintf('%-10s','c (uM^-1)'); fprintf('%12s',pars{:}); fprintf('\n');
for i=1:numel(cv)
  fprintf('%-10.4f',cv(i)); fprintf('%12.4f',Sens(i,:)); fprintf('\n');
end

% rank correlation of the sensitivity ordering across calibrations
[~,r1]=sort(abs(Sens(1,:)),'descend');
[~,r2]=sort(abs(Sens(end,:)),'descend');
fprintf('\nparameter ranking, c = %.4f : %s\n',cv(1),strjoin(pars(r1),' > '));
fprintf('parameter ranking, c = %.4f : %s\n',cv(end),strjoin(pars(r2),' > '));

figure(1)
plot(tterm/yr,cv,'-rs','LineWidth',1.5,'MarkerSize',9,...
    'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
hold on
plot(6,c17,'kp','MarkerSize',16,'MarkerFaceColor','k');
hold off
xlabel('assumed survival after diagnosis used for calibration (years)');
ylabel('implied c (\muM^{-1})');
legend('implied by Eq. (51)','value used in this study','location','best');

figure(2)
plot(cv,abs(Sens(:,1)),'-rs','LineWidth',1.5);
hold on
plot(cv,abs(Sens(:,2)),'--bdiamond','LineWidth',1.5);
plot(cv,abs(Sens(:,4)),'-.g^','LineWidth',1.5);
plot(cv,abs(Sens(:,7)),':mv','LineWidth',1.5);
hold off
xlabel('c (\muM^{-1})');
ylabel('|dimensionless sensitivity| at 6 years');
legend('[H]_{tot}','k_{cl}','\theta_{1/2,B}','k_1','location','best');

set(figure(1),'Name','Fig. S28a','NumberTitle','off');
savefig_emf(figure(1),'figS28a');
set(figure(2),'Name','Fig. S28b','NumberTitle','off');
savefig_emf(figure(2),'figS28b');

toc
