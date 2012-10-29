function CIeng2PT
%This is a generic eng to PT file conversion program.
%A specific program for each instrument/axis will need to be written.
%Probably best to put it into the axis directory.
% CIeng2PT
% Reads tildeeng_1.mat and creates a HHH-specific PT.mat appropriate
% for use with the icosfit utilities.
E = load('tildeeng_1.mat');
E16 = load('tildeeng_16.mat');
PT.TPT = E.Ttildeeng_1;
PT.ScanNum = E.SSP3_Num;
PT.QCLI_Wave = E.QCLI3_Wave;
PT.CellP = E16.HCInP - 0.0;
PT.Tavg = 273.15 + E.HCExT;
save PT.mat -STRUCT PT
