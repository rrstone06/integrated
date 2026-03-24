clear; clc; close all;
%% Format Data

data = readtable('integrated.csv');
data = data(:, 4:6);
data.Properties.VariableNames = {'time', 'value', 'name'};

LC1 = data(strcmp(data.name,'LC_1'), {'time','value'}); LC1.Properties.VariableNames{2} = 'LC1';
LC2 = data(strcmp(data.name,'LC_2'), {'time','value'}); LC2.Properties.VariableNames{2} = 'LC2';
LC3 = data(strcmp(data.name,'LC_3'), {'time','value'}); LC3.Properties.VariableNames{2} = 'LC3';
LC4 = data(strcmp(data.name,'LC_4'), {'time','value'}); LC4.Properties.VariableNames{2} = 'LC4';
PT  = data(strcmp(data.name,'PT_P_E_F_7'), {'time','value'}); PT.Properties.VariableNames{2} = 'PT';

% Join all on timestamp
joined = innerjoin(LC1, LC2, 'Keys', 'time');
joined = innerjoin(joined, LC3, 'Keys', 'time');
joined = innerjoin(joined, LC4, 'Keys', 'time');
joined = innerjoin(joined, PT,  'Keys', 'time');

% Remove rows where any value is zero
mask = ~any(joined{:, 2:end} == 0, 2);
joined = joined(mask, :);

% Final matrices
LC = joined{:, {'LC1','LC2','LC3','LC4'}};
PT = joined.PT;