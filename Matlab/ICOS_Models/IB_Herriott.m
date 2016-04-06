%%
% IB Herriott Power investigation
% load an IB file, then let's try to come up with better metrics on the
% performance.
% [ ]Expected power can be based on advancement angle, 180/asind(Rw1/Rr1)
%   Except, duh, we don't know Rw1 and Rr1
% [ ]Calculate losses on first bounce (i.e. out the aperature)
% [ ]Calculate losses outside on RIM
% [ ]Calculate losses outside on M1
datadir = 'C:\Data\ES96\ICOS\3D\HCl_TB4';
load([datadir filesep 'IB_HCl_TB4_rd08_th13.0.393_75_-281_71_26.9.19_50x100.mat']);
P = IB.P;
IB2 = ICOS_beam(@ICOS_Model6, P);
IB2.Sample('beam_samples', 100, ...
  'ICOS_passes', 50, 'opt_n', 6, ...
  'n_optics', 6, 'Track_Power', 1, ...
  'Herriott_passes', 100, ...
  'mnc', 'test');
%%
R0_2 = IB2.Res.Pwr.R0_2;
R2_1 = IB2.Res.Pwr.R2_1;
R1_2 = IB2.Res.Pwr.R1_2;
NS = R0_2.NI + R0_2.NO;
PI = (R0_2.I + R0_2.O)/NS; % should be 1
PO_E1 = R2_1.E1/NS;
PO_EO = R2_1.EO/NS;
PO_1O = R2_1.O/NS;
PO_2O = (R0_2.O + R1_2.O)/NS;
PO_ICOS = (R0_2.I + R1_2.I)*IB2.P.T/NS;
PO_Hloss = R2_1.I * (1-IB2.P.HR)/NS;
PO = PO_E1+PO_EO+PO_1O+PO_2O+PO_Hloss+PO_ICOS;
NH = (R0_2.I+R1_2.I)/(R0_2.I+R0_2.O);
%%
% Average number of passes for spots that exit through aperature after
% more than one bounce:
N = log(PO_EO*NS/(R2_1.NEO*PI*(1-IB2.P.T)))/log((1-IB2.P.T)*IB2.P.HR) + 1;
%%
% Now if every beam resulted in N passes thorugh Herriott cell, we'd
% get N spots on the back of the ICOS mirror, and a total power on the
% back of the ICOS mirror of:
r = (1-IB2.P.T)*IB2.P.HR;
P2_ideal = (1-r^N)/(1-r);
%%
dNP = diff(IB2.Res.NPass(1:IB2.Res.N+1,:));
v = find(dNP(:,1)<0 | (dNP(:,1)==0 & dNP(:,2)<0));
NPH = NP(v,1);
NH = mode(NPH);
