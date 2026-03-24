
%Assumes Data for Each load cell is organized by column (i.e. all data for
%load cell 1 is column 1, all data for load cell 2 is column 2, etc.)

%Assumes all of this data is conbined into an X by 5 matrix, where each
%column corresponds to its respective load cell(with the 5th corresponding
%to the total force value readout, and X represents the
%number of data entries each load cell has. Also assumes every entry in
%this matrix, stored as cellData, has no NaN or 0 data values.

clear; clc; close all;
%% Format Data
data = readtable('integrated.csv');
data = data(:, 4:6);
data.Properties.VariableNames = {'time', 'value', 'name'};
LC1 = data(strcmp(data.name,'LC_1'), {'time','value'}); LC1.Properties.VariableNames{2} = 'LC1';
LC2 = data(strcmp(data.name,'LC_2'), {'time','value'}); LC2.Properties.VariableNames{2} = 'LC2';
LC3 = data(strcmp(data.name,'LC_3'), {'time','value'}); LC3.Properties.VariableNames{2} = 'LC3';
LC4 = data(strcmp(data.name,'LC_4'), {'time','value'}); LC4.Properties.VariableNames{2} = 'LC4';
PT = data(strcmp(data.name,'PT_P_E_F_7'), {'time','value'}); PT.Properties.VariableNames{2} = 'PT';
% Join all on timestamp
joined = innerjoin(LC1, LC2, 'Keys', 'time');
joined = innerjoin(joined, LC3, 'Keys', 'time');
joined = innerjoin(joined, LC4, 'Keys', 'time');
joined = innerjoin(joined, PT, 'Keys', 'time');
% Remove rows where any value is zero
mask = ~any(joined{:, 2:end} == 0, 2);
joined = joined(mask, :);
% Final matrices
cellData = joined{:, {'LC1','LC2','LC3','LC4'}};
PT = joined.PT.*4.98;


lc1Coeff = []; %matrix storing coeffecients for load cells 1-4
lc2Coeff = [];
lc3Coeff = [];
lc4Coeff = [];
i = 1;
while i < length(cellData)
    
eq1 = [cellData(i,1) cellData(i,2) cellData(i,3) cellData(i,4)];
eq2 = [cellData(i+1,1) cellData(i+1,2) cellData(i+1,3) cellData(i+1,4)];
eq3 = [cellData(i+2,1) cellData(i+2,2) cellData(i+2,3) cellData(i+2,4)];
eq4 = [cellData(i+3,1) cellData(i+3,2) cellData(i+3,3) cellData(i+3,4)];

forceTotals = [PT(i); PT(i+1); PT(i+2); PT(i+3)];

eqns = [eq1;eq2;eq3;eq4];

tempCoeff = inv(eqns)*forceTotals;
lc1Coeff = [lc1Coeff; tempCoeff(1)];
lc2Coeff = [lc2Coeff; tempCoeff(2)];
lc3Coeff = [lc3Coeff; tempCoeff(3)];
lc4Coeff = [lc4Coeff; tempCoeff(4)];

i = i+4;

end
xAxis = (1:1:507);
figure(1)
plot(xAxis,lc1Coeff)
hold on
plot(xAxis,lc2Coeff)
plot(xAxis,lc3Coeff)
plot(xAxis,lc4Coeff)