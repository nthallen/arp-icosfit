function CIeng2PT
% HHHeng2PT
% Reads tildeeng_1.mat and creates a HHH-specific PT.mat appropriate
% for use with the icosfit utilities.
E = load('tildeeng_1.mat');
E16 = load('tildeeng_16.mat');
PT.TPT = E.Ttildeeng_1;
PT.CPCI14 = E.SSP3_Num;
PT.QCLI_Wave = E.QCLI3_Wave;
PT.CellP = E16.HCInP - 0.0;
PT.Tavg = 273.15 + E.HCExT;
save PT.mat -STRUCT PT
