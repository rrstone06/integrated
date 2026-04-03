%Assumes Data for Each load cell is organized by column (i.e. all data for
%load cell 1 is column 1, all data for load cell 2 is column 2, etc.)
%Assumes all of this data is combined into an X by 5 matrix, where each
%column corresponds to its respective load cell (with the 5th corresponding
%to the total force value readout, and X represents the
%number of data entries each load cell has. Also assumes every entry in
%this matrix, stored as cellData, has no NaN or 0 data values.
clear; clc; close all;

%% Format Data
data = readtable('integrated_1ms.csv');
data = data(:, [6:7, 13]);
data.Properties.VariableNames = {'time', 'value', 'name'};

LC1 = data(strcmp(data.name,'LC_1'), {'time','value'}); LC1.Properties.VariableNames{2} = 'LC1';
LC2 = data(strcmp(data.name,'LC_2'), {'time','value'}); LC2.Properties.VariableNames{2} = 'LC2';
LC3 = data(strcmp(data.name,'LC_3'), {'time','value'}); LC3.Properties.VariableNames{2} = 'LC3';
LC4 = data(strcmp(data.name,'LC_4'), {'time','value'}); LC4.Properties.VariableNames{2} = 'LC4';
PT  = data(strcmp(data.name,'PT_P_E_F_7'), {'time','value'}); PT.Properties.VariableNames{2} = 'PT';

% Convert timestamps to datetime
parseDT = @(c) datetime(regexprep(string(c), '(\d{2}:\d{2}:\d{2})Z', '$1.000Z'), ...
               'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSSZ', 'TimeZone', 'UTC');

t1  = parseDT(LC1.time);
t2  = parseDT(LC2.time);
t3  = parseDT(LC3.time);
t4  = parseDT(LC4.time);
tPT = parseDT(PT.time);

% Convert to seconds
t1s  = seconds(t1  - t1(1));
t2s  = seconds(t2  - t1(1));
t3s  = seconds(t3  - t1(1));
t4s  = seconds(t4  - t1(1));
tPTs = seconds(tPT - t1(1));

% Match all sensors to LC1 timestamps within 10ms tolerance
tol = 0.01;
keepIdx = false(length(t1s), 1);
lc2Val = zeros(length(t1s),1); lc3Val = zeros(length(t1s),1);
lc4Val = zeros(length(t1s),1); ptVal  = zeros(length(t1s),1);

for i = 1:length(t1s)
    d2 = abs(t2s - t1s(i)); [m2, i2] = min(d2);
    d3 = abs(t3s - t1s(i)); [m3, i3] = min(d3);
    d4 = abs(t4s - t1s(i)); [m4, i4] = min(d4);
    dP = abs(tPTs - t1s(i)); [mP, iP] = min(dP);

    if m2 <= tol && m3 <= tol && m4 <= tol && mP <= tol
        keepIdx(i)  = true;
        lc2Val(i)   = LC2.LC2(i2);
        lc3Val(i)   = LC3.LC3(i3);
        lc4Val(i)   = LC4.LC4(i4);
        ptVal(i)    = PT.PT(iP);
    end
end

% Trim to matched rows only
tGrid  = t1s(keepIdx);
lc1Val = LC1.LC1(keepIdx);
lc2Val = lc2Val(keepIdx);
lc3Val = lc3Val(keepIdx);
lc4Val = lc4Val(keepIdx);
ptVal  = ptVal(keepIdx);

cellData = [lc1Val, lc2Val, lc3Val, lc4Val];
PT_vec   = ptVal .* 4.98;

fprintf('Matched rows: %d\n', sum(keepIdx));

% Remove rows where any LC or PT is zero
mask     = ~any(cellData == 0, 2) & (PT_vec ~= 0);
cellData = cellData(mask, :);
PT_vec   = PT_vec(mask);
tGrid    = tGrid(mask);

%% Solve as Overdetermined Least Squares (all data at once)
coeffs = cellData \ PT_vec;

fprintf('LC1 coefficient: %.6f\n', coeffs(1));
fprintf('LC2 coefficient: %.6f\n', coeffs(2));
fprintf('LC3 coefficient: %.6f\n', coeffs(3));
fprintf('LC4 coefficient: %.6f\n', coeffs(4));

%% Least Squares by Thrust Range
binEdges  = linspace(0, 10500, 8);
binLabels = {};
rangeCoeffs = zeros(7, 4);

for b = 1:7
    lo = binEdges(b);
    hi = binEdges(b+1);

    binMask = PT_vec >= lo & PT_vec < hi;
    binLabels{b} = sprintf('%.0f-%.0f lbf', lo, hi);

    if sum(binMask) < 4
        fprintf('Bin %d (%.0f-%.0f lbf): not enough data\n', b, lo, hi);
        continue
    end

    binCoeffs = cellData(binMask, :) \ PT_vec(binMask);
    rangeCoeffs(b, :) = binCoeffs';

    fprintf('Bin %d (%.0f-%.0f lbf):\n', b, lo, hi);
    fprintf('  LC1: %.6f  LC2: %.6f  LC3: %.6f  LC4: %.6f\n', ...
        binCoeffs(1), binCoeffs(2), binCoeffs(3), binCoeffs(4));
end

%% Plot LC values and PT over time
figure(1)
plot(tGrid, cellData(:,1))
hold on
plot(tGrid, cellData(:,2))
plot(tGrid, cellData(:,3))
plot(tGrid, cellData(:,4))
plot(tGrid, PT_vec, 'k--')
legend('LC1 (lbf)','LC2 (lbf)','LC3 (lbf)','LC4 (lbf)','PT (lbf)')
xlabel('Time (s)')
ylabel('Force (lbf)')
title('Raw LC and PT Values Over Time')

%% Plot Coefficients by Thrust Range
lcNames = {'LC1','LC2','LC3','LC4'};
figure(2)
tiledlayout(4, 1)
for lc = 1:4
    nexttile
    bar(rangeCoeffs(:, lc))
    xticklabels(binLabels)
    xtickangle(30)
    ylabel('Coefficient (dimensionless)')
    title(sprintf('%s Coefficient by Thrust Range', lcNames{lc}))
end
sgtitle('Load Cell Coefficients by Thrust Range (lbf)')

%% Test Coefficients
testCoeffs(coeffs', 'integrated_1ms.csv')