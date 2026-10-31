clc, clear, close all

% Reproduces Fig. 1 of the paper: (a) fraction of lipid-free SAA, omegaS, and
% (b) fraction of SAA-binding sites on HDL that are unoccupied, omegaH, as
% functions of [H]tot, from the SAA-HDL binding equilibrium (Section 2.1.1).
% Two cases: the baseline Kd and the tight-binding limit Kd = 1e-10 uM.
% Standalone; no ODEs are solved.

Htot=linspace(1,30,100);
Stot=8;

% Baseline Kd: the value at which omegaS = 0.05 for [H]tot = 30 uM and
% Stot = 8 uM, i.e. Kd = omegaS*(Htot-Stot+omegaS*Stot)/(1-omegaS).
Kd=1.1789;

Sfree=0.5*(-(Kd+Htot-Stot)+sqrt((Kd+Htot-Stot).^2+4*Kd*Stot));
omegaS=Sfree./Stot;
Hfree=Htot-Stot+Sfree;
omegaH=Hfree./Htot;

omegaS_baselineKd=omegaS;
omegaH_baselineKd=omegaH;

Kd=1e-10;

Sfree=0.5*(-(Kd+Htot-Stot)+sqrt((Kd+Htot-Stot).^2+4*Kd*Stot));
omegaS=Sfree./Stot;
Hfree=Htot-Stot+Sfree;
omegaH=Hfree./Htot;

omegaS_smallKd=omegaS;
omegaH_smallKd=omegaH;

figure('Name','Fig. 1a','NumberTitle','off')
plot(Htot,omegaS_baselineKd,'-r','LineWidth',1.5)
hold on
plot(Htot,omegaS_smallKd,'--b','LineWidth',1.5)
hold off
xlabel('[H]_{tot}  (\muM)');
ylabel('\omega_S');
legend('K_d=1.1789 \muM','K_d=10^{-10} \muM','location','best')

figure('Name','Fig. 1b','NumberTitle','off')
plot(Htot,omegaH_baselineKd,'-r','LineWidth',1.5)
hold on
plot(Htot,omegaH_smallKd,'--b','LineWidth',1.5)
hold off
xlabel('[H]_{tot}  (\muM)');
ylabel('\omega_H');
legend('K_d=1.1789 \muM','K_d=10^{-10} \muM','location','best')
