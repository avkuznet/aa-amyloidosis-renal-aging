% Baseline parameter values of Table 2, shared by AA_tipping_point.m,
% AA_flares.m, AA_treatment.m, AA_validation.m and AA_calibration_c.m.
% The exception is Htot: these scripts use [H]tot = 15 uM (Figs. 8-11, S27
% and the calibration of c), whereas the Table 3 baseline is 30 uM, which
% AA_calibration_c.m sets explicitly for its sensitivity calculation.
% The parameter-study scripts (AA_Htot.m, AA_kcl.m, ...) set the same
% values explicitly.

Vkidneys=4.04e-4;
Akidneys=0.6;
hS=5.55e-7;
kcl=9.65e-5;
k1=1e-4;
k2=1e-4;
theta12B=8.64e4;

Kd=1.1789;
Stot=8;
Htot=15;

MWAA=1.2e4;
a21=1e21;
rofibril=1.35e-12;

omega=1.4e10;
c17=0.334;
tAAonset=50*(365*24*3600);
yr=365*24*3600;

% functional protein clearance
T12A_f=3.6e3; T12B_f=1.73e5; T12I_f=3.16e7; T12S_f=3.6e3;
% impaired protein degradation pathways
T12A_i=1e14;  T12B_i=1e14;  T12I_i=1e14;  T12S_i=1e14;
