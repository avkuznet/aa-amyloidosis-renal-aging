clc, clear, close all
tic
% Reproduces Fig. 8 of the paper (Section 2.5).
% Autocatalytic tipping point of the F-W step of the AA model.
% Reduces the steady-state system (21)-(22) to a single dimensionless
% relation b = B/B* governed by kappa1 and kappa2, and locates the
% transcritical threshold kappa2 = 1.

AA_params

T12A=T12A_f; T12B=T12B_f; T12I=T12I_f; T12S=T12S_f;

bfun = @(k1d,k2d) (-(1+k1d-k2d)+sqrt((1+k1d-k2d).^2+4*k1d.*k2d))./(2*k2d);

%%%%%%%%%%%%%%%% panel (a): b vs kappa2 for three kappa1
kap2v=linspace(0.01,4,400);
kap1v=[1e-3 1e-2 1e-1];

figure(1)
plot(kap2v,bfun(kap1v(1),kap2v),'-r','LineWidth',1.5);
hold on
plot(kap2v,bfun(kap1v(2),kap2v),'--b','LineWidth',1.5);
plot(kap2v,bfun(kap1v(3),kap2v),'-.g','LineWidth',1.5);
plot(kap2v,max(0,(kap2v-1)./kap2v),':k','LineWidth',2);
plot([1 1],[0 1],':m','LineWidth',1.5);
hold off
xlabel('\kappa_2');
ylabel('b = [B]/[B]^*');
legend('\kappa_1=10^{-3}','\kappa_1=10^{-2}','\kappa_1=10^{-1}',...
    '\kappa_1\rightarrow0 (1-1/\kappa_2)','\kappa_2=1','location','best');
ylim([0 1])

%%%%%%%%%%%%%%%% panel (b): b vs Htot for three k2, with ODE verification
Htotv=linspace(5,40,200);
k2v=[1e-6 1e-5 1e-4];
bH=zeros(numel(k2v),numel(Htotv));
kap2H=zeros(numel(k2v),numel(Htotv));
Hcrit=nan(1,numel(k2v));

for j=1:numel(k2v)
  for i=1:numel(Htotv)
    H=Htotv(i);
    Sfree=0.5*(-(Kd+H-Stot)+sqrt((Kd+H-Stot).^2+4*Kd*Stot));
    Sss=Akidneys*hS*Sfree*T12S/(Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));
    Bstar=kcl*Sss*T12B*theta12B/(log(2)*(T12B+theta12B));
    kap1=k1*T12A/log(2);
    kap2=k2v(j)*Bstar*T12A/log(2);
    kap2H(j,i)=kap2;
    bH(j,i)=bfun(kap1,kap2);
  end
  % critical HDL concentration at which kappa2 = 1
  f=@(H) k2v(j)*(kcl*(Akidneys*hS* ...
      (0.5*(-(Kd+H-Stot)+sqrt((Kd+H-Stot).^2+4*Kd*Stot))) *T12S/ ...
      (Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2))) ...
      *T12B*theta12B/(log(2)*(T12B+theta12B)))*T12A/log(2)-1;
  try
    Hcrit(j)=fzero(f,[1e-3 500]);
  catch
    Hcrit(j)=NaN;
  end
  fprintf('k2=%8.1e  kappa2(Htot=15)=%9.4f  kappa2(Htot=30)=%9.4f  Hcrit=%8.3f uM\n', ...
      k2v(j),kap2H(j,Htotv==15 | abs(Htotv-15)==min(abs(Htotv-15))),...
      kap2H(j,abs(Htotv-30)==min(abs(Htotv-30))),Hcrit(j));
end

% verify the reduced solution against ode45 at a few HDL values
Hchk=[8 15 22 30 38];
Bode=zeros(size(Hchk)); Bred=zeros(size(Hchk)); bode=zeros(size(Hchk));
options=odeset('RelTol',1e-12,'AbsTol',1e-14);
for i=1:numel(Hchk)
  H=Hchk(i);
  Sfree=0.5*(-(Kd+H-Stot)+sqrt((Kd+H-Stot).^2+4*Kd*Stot));
  sol=ode45(@(t,y) funCDCMCB_AA(t,y,hS,Akidneys,Vkidneys,Sfree,kcl,T12S, ...
      k1,k2,T12A,T12B,theta12B,T12I),[0 1e7],[0 0 0 0],options);
  y=deval(sol,1e7);
  Sss=Akidneys*hS*Sfree*T12S/(Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));
  Bstar=kcl*Sss*T12B*theta12B/(log(2)*(T12B+theta12B));
  kap1=k1*T12A/log(2); kap2=k2*Bstar*T12A/log(2);
  Bode(i)=y(3); Bred(i)=bfun(kap1,kap2)*Bstar; bode(i)=y(3)/Bstar;
  fprintf('Htot=%5.1f  B_ode=%12.6e  B_reduced=%12.6e  rel.err=%9.2e\n',...
      H,Bode(i),Bred(i),abs(Bode(i)-Bred(i))/Bode(i));
end

figure(2)
semilogy(Htotv,bH(1,:),'-r','LineWidth',1.5);
hold on
semilogy(Htotv,bH(2,:),'--b','LineWidth',1.5);
semilogy(Htotv,bH(3,:),'-.g','LineWidth',1.5);
semilogy(Hchk,bode,'ks','LineWidth',1.5,'MarkerSize',9,'MarkerFaceColor',[.8 .8 .8]);
hold off
xlabel('[H]_{tot}  (\muM)');
ylabel('b = [B]/[B]^*');
legend('k_2=10^{-6} \muM^{-1} s^{-1}','k_2=10^{-5} \muM^{-1} s^{-1}',...
    'k_2=10^{-4} \muM^{-1} s^{-1}','numerical solution of Eqs. (15),(19),(20)',...
    'location','best');

fprintf('\nkappa1 (baseline k1=1e-4, functional clearance) = %8.4f\n',k1*T12A/log(2));

%%%%%%%%%%%%%%%% linear stability of the monomer-oligomer subsystem
% Eqs. (19)-(20) at fixed renal SAA:  dA/dt = kcl*S - k1*A - k2*A*B - lamA*A
%                                    dB/dt = k1*A + k2*A*B - lamB*B
% For kappa1 = 0 the zero-oligomer state has eigenvalues -lamA and
% lamB*(kappa2-1), so it loses stability at kappa2 = 1, where it exchanges
% stability with the branch b = 1-1/kappa2; that branch exists for every
% kappa2 but is negative (unphysical) below the threshold and stable above
% it (transcritical bifurcation).
% For kappa1 > 0 the single non-negative steady state is stable throughout.
lamA=log(2)/T12A;
lamB=log(2)/T12B+log(2)/theta12B;
H=15;
Sfree=0.5*(-(Kd+H-Stot)+sqrt((Kd+H-Stot).^2+4*Kd*Stot));
Sss=Akidneys*hS*Sfree*T12S/(Akidneys*hS*T12S+kcl*T12S*Vkidneys+Vkidneys*log(2));
Bstar=kcl*Sss/lamB;
kap2s=[0.25 0.5 0.9 1.1 2 4];
fprintf('\nlinear stability, kappa1 = 0 (k1 = 0)\n');
fprintf(' kappa2   maxRe(eig) zero branch   lamB*(kappa2-1)   maxRe(eig) positive branch\n');
for kap2=kap2s
  k2s=kap2*lamA/Bstar;                      % k2 that gives this kappa2
  A0=kcl*Sss/lamA;                          % zero-oligomer state
  J0=[-lamA, -k2s*A0; 0, k2s*A0-lamB];
  e0=max(real(eig(J0)));
  if kap2>1
    A1=lamB/k2s; B1=(kcl*Sss-lamA*A1)/(k2s*A1);
    J1=[-k2s*B1-lamA, -k2s*A1; k2s*B1, k2s*A1-lamB];
    fprintf('%7.2f  %+24.4e  %+16.4e  %+24.4e\n',kap2,e0,lamB*(kap2-1),max(real(eig(J1))));
  else
    fprintf('%7.2f  %+24.4e  %+16.4e  %24s\n',kap2,e0,lamB*(kap2-1),'(b < 0, unphysical)');
  end
end
fprintf('\nlinear stability, baseline kappa1 = %6.4f\n',k1*T12A/log(2));
fprintf(' kappa2        b    maxRe(eig)\n');
for kap2=kap2s
  k2s=kap2*lamA/Bstar;
  bb=bfun(k1/lamA,kap2); B1=bb*Bstar;
  A1=kcl*Sss/(k1+k2s*B1+lamA);
  J=[-k1-k2s*B1-lamA, -k2s*A1; k1+k2s*B1, k2s*A1-lamB];
  fprintf('%7.2f  %7.4f  %+12.4e\n',kap2,bb,max(real(eig(J))));
end

set(figure(1),'Name','Fig. 8a','NumberTitle','off');
savefig_emf(figure(1),'fig8a');
set(figure(2),'Name','Fig. 8b','NumberTitle','off');
savefig_emf(figure(2),'fig8b');

toc
