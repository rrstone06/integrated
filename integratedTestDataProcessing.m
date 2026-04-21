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

% Use LC1 as reference, find indices where all others have a real point within tolerance
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
PT_vec   = ptVal .* (((2.52/2)^2)*pi);

fprintf('Matched rows: %d\n', sum(keepIdx));

% Remove rows where any LC or PT is zero
mask     = ~any(cellData == 0, 2) & (PT_vec ~= 0);
cellData = cellData(mask, :);
PT_vec   = PT_vec(mask);
tGrid    = tGrid(mask);

% Time mask — exclude idle period between tests
tMin = 1200;  % seconds
tMax = 1400;  % seconds
timeMask = ~(tGrid >= tMin & tGrid <= tMax);
cellData = cellData(timeMask, :);
PT_vec   = PT_vec(timeMask);
tGrid    = tGrid(timeMask);
fprintf('Rows after time mask: %d\n', sum(timeMask));

%% Split into peaks by time gaps
dt        = diff(tGrid);
gapThresh = 1.0;
gapIdx    = find(dt > gapThresh);

segStarts = [1; gapIdx+1];
segEnds   = [gapIdx; length(tGrid)];

fprintf('Found %d peaks\n', length(segStarts));

peakCoeffs = zeros(length(segStarts), 5);

for p = 1:length(segStarts)
    idx   = segStarts(p):segEnds(p);
    cdSeg = cellData(idx, :);
    ptSeg = PT_vec(idx);
    tSeg  = tGrid(idx);

    fprintf('\n--- Peak %d (%d rows, t=%.1f to %.1f s) ---\n', ...
        p, length(idx), tSeg(1), tSeg(end));

    % Initial solve
    cb = [cdSeg, ones(size(cdSeg,1),1)] \ ptSeg;

    % Remove outliers
    res   = ptSeg - (cdSeg * cb(1:4) + cb(5));
    omask = abs(res - mean(res)) <= 3 * std(res);
    cdSeg = cdSeg(omask, :);
    ptSeg = ptSeg(omask);

    fprintf('Removed %d outliers (%.1f%%)\n', sum(~omask), 100*mean(~omask));

    % Re-solve
    cb = [cdSeg, ones(size(cdSeg,1),1)] \ ptSeg;
    peakCoeffs(p,:) = cb';

    % Error metrics
    est = cdSeg * cb(1:4) + cb(5);
    fprintf('LC1: %.6f  LC2: %.6f  LC3: %.6f  LC4: %.6f  Bias: %.4f\n', ...
        cb(1), cb(2), cb(3), cb(4), cb(5));
    fprintf('RMSE: %.4f lbf  Max: %.4f lbf  Mean: %.4f lbf\n', ...
        sqrt(mean((ptSeg-est).^2)), max(abs(ptSeg-est)), mean(abs(ptSeg-est)));

    % Plot peak fit
    figure(2 + p)
    plot(tSeg(omask), ptSeg, 'k', 'DisplayName', 'PT (actual)')
    hold on
    plot(tSeg(omask), est, 'r--', 'DisplayName', 'LC estimate')
    plot(tSeg(omask), ptSeg - est, 'b', 'DisplayName', 'Residual')
    yline(0, 'k--')
    legend
    xlabel('Time (s)')
    ylabel('Force (lbf)')
    title(sprintf('Peak %d Fit  |  RMSE: %.1f lbf  Mean: %.1f lbf', ...
        p, sqrt(mean((ptSeg-est).^2)), mean(abs(ptSeg-est))))
end

%% Average coefficients across peaks
if size(peakCoeffs, 1) == 1
    coeffs = peakCoeffs(1, 1:4)';
    bias   = peakCoeffs(1, 5);
else
    coeffs = mean(peakCoeffs(:,1:4))';
    bias   = mean(peakCoeffs(:,5));
end
fprintf('\n--- Average across peaks ---\n');
fprintf('LC1: %.6f  LC2: %.6f  LC3: %.6f  LC4: %.6f  Bias: %.4f lbf\n', ...
    coeffs(1), coeffs(2), coeffs(3), coeffs(4), bias);

%% Solve by Thrust Range
binEdges    = linspace(5500, 10500, 6);
binLabels   = {};
rangeCoeffs = zeros(5, 5);

for b = 1:5
    lo = binEdges(b);
    hi = binEdges(b+1);

    binMask = PT_vec >= lo & PT_vec < hi;
    binLabels{b} = sprintf('%.0f-%.0f lbf', lo, hi);

    if sum(binMask) < 4
        fprintf('Bin %d (%.0f-%.0f lbf): not enough data\n', b, lo, hi);
        continue
    end

    binData   = cellData(binMask, :);
    binCoeffs = [binData, ones(size(binData,1),1)] \ PT_vec(binMask);
    rangeCoeffs(b, :) = binCoeffs';

    fprintf('Bin %d (%.0f-%.0f lbf):\n', b, lo, hi);
    fprintf('  LC1: %.6f  LC2: %.6f  LC3: %.6f  LC4: %.6f  Bias: %.4f\n', ...
        binCoeffs(1), binCoeffs(2), binCoeffs(3), binCoeffs(4), binCoeffs(5));
end

%% Plot LC values and PT over time
dt_plot = diff(tGrid);
gapMask = [false; dt_plot > 1.0];

tPlot  = tGrid;
cdPlot = cellData;
ptPlot = PT_vec;

tPlot(gapMask)    = NaN;
cdPlot(gapMask,:) = NaN;
ptPlot(gapMask)   = NaN;

figure(1)
plot(tPlot, cdPlot(:,1))
hold on
plot(tPlot, cdPlot(:,2))
plot(tPlot, cdPlot(:,3))
plot(tPlot, cdPlot(:,4))
plot(tPlot, ptPlot, 'k--')
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
testCoeffs([coeffs; bias]', 'integrated_1ms.csv')