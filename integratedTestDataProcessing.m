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

% Convert to seconds (numeric) for interpolation
t1s  = seconds(t1  - t1(1));
t2s  = seconds(t2  - t1(1));
t3s  = seconds(t3  - t1(1));
t4s  = seconds(t4  - t1(1));
tPTs = seconds(tPT - t1(1));

% Use LC1 as reference, find indices where all others have a real point within tolerance
tol = 0.01; % 10ms tolerance
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
t1s    = t1s(keepIdx);
lc1Val = LC1.LC1(keepIdx);
lc2Val = lc2Val(keepIdx);
lc3Val = lc3Val(keepIdx);
lc4Val = lc4Val(keepIdx);
ptVal  = ptVal(keepIdx);

cellData = [lc1Val, lc2Val, lc3Val, lc4Val];
PT_vec   = ptVal .* 4.98;
tGrid    = t1s;

fprintf('Matched rows: %d\n', sum(keepIdx));

% Remove rows where any LC or PT is zero
mask     = ~any(cellData == 0, 2) & (PT_vec ~= 0);
cellData = cellData(mask, :);
PT_vec   = PT_vec(mask);
tGrid    = tGrid(mask);

%% Solve for Coefficients
lc1Coeff = [];
lc2Coeff = [];
lc3Coeff = [];
lc4Coeff = [];
skipped  = 0;

step      = 1;   % rows apart within a group
groupStep = 4;   % move to next group every 4 rows
i = 1;
while i + 3*step <= length(cellData)
    idx = [i, i+step, i+2*step, i+3*step];
    eqns        = cellData(idx, :);
    forceTotals = PT_vec(idx);

    if rcond(eqns) < 1e-12
        skipped = skipped + 1;
        i = i + groupStep;
        continue
    end

    tempCoeff = inv(eqns) * forceTotals;

    lc1Coeff = [lc1Coeff; tempCoeff(1)];
    lc2Coeff = [lc2Coeff; tempCoeff(2)];
    lc3Coeff = [lc3Coeff; tempCoeff(3)];
    lc4Coeff = [lc4Coeff; tempCoeff(4)];

    i = i + groupStep;
end

fprintf('Skipped %d singular groups out of %d total\n', skipped, floor(length(cellData)/groupStep));

%% Floating Average for Each Coefficient
windowSize = 50;

lc1Avg = movmean(lc1Coeff, windowSize);
lc2Avg = movmean(lc2Coeff, windowSize);
lc3Avg = movmean(lc3Coeff, windowSize);
lc4Avg = movmean(lc4Coeff, windowSize);

fprintf('LC1 mean coefficient: %.6f\n', mean(lc1Coeff));
fprintf('LC2 mean coefficient: %.6f\n', mean(lc2Coeff));
fprintf('LC3 mean coefficient: %.6f\n', mean(lc3Coeff));
fprintf('LC4 mean coefficient: %.6f\n', mean(lc4Coeff));

%% Floating Average by Thrust Range
binEdges = linspace(0, 10500, 8);
binLabels = {};
rangeCoeffs = zeros(7, 4);

ptGroups = PT_vec(1:groupStep:end);
ptGroups = ptGroups(1:length(lc1Coeff));
x = (1:length(lc1Coeff))';

for b = 1:7
    lo = binEdges(b);
    hi = binEdges(b+1);

    binMask = ptGroups >= lo & ptGroups < hi;
    binLabels{b} = sprintf('%.0f-%.0f lbf', lo, hi);

    if sum(binMask) < 2
        fprintf('Bin %d (%.0f-%.0f lbf): not enough data\n', b, lo, hi);
        continue
    end

    rangeCoeffs(b, :) = [mean(lc1Coeff(binMask)), mean(lc2Coeff(binMask)), ...
                          mean(lc3Coeff(binMask)), mean(lc4Coeff(binMask))];

    fprintf('Bin %d (%.0f-%.0f lbf):\n', b, lo, hi);
    fprintf('  LC1: %.6f  LC2: %.6f  LC3: %.6f  LC4: %.6f\n', ...
        rangeCoeffs(b,1), rangeCoeffs(b,2), rangeCoeffs(b,3), rangeCoeffs(b,4));
end

%% Plot Coefficients Over Time
tCoeff = tGrid(1:groupStep:end);
tCoeff = tCoeff(1:length(lc1Coeff));

figure(1)
plot(tCoeff, lc1Coeff)
hold on
plot(tCoeff, lc2Coeff)
plot(tCoeff, lc3Coeff)
plot(tCoeff, lc4Coeff)
plot(tCoeff, lc1Avg, 'b--')
plot(tCoeff, lc2Avg, 'r--')
plot(tCoeff, lc3Avg, 'g--')
plot(tCoeff, lc4Avg, 'm--')
legend('LC1','LC2','LC3','LC4','LC1 avg','LC2 avg','LC3 avg','LC4 avg')
xlabel('Time (s)')
ylabel('Coefficient (dimensionless)')
title('Load Cell Coefficients Over Time')

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
bestCoeffs = [mean(lc1Coeff), mean(lc2Coeff), mean(lc3Coeff), mean(lc4Coeff)];
testCoeffs(bestCoeffs, 'integrated_1ms.csv')