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

parseDT = @(c) datetime(regexprep(string(c), '(\d{2}:\d{2}:\d{2})Z', '$1.000Z'), ...
               'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSSZ', 'TimeZone', 'UTC');

t1  = parseDT(LC1.time);
t2  = parseDT(LC2.time);
t3  = parseDT(LC3.time);
t4  = parseDT(LC4.time);
tPT = parseDT(PT.time);

t1s  = seconds(t1  - t1(1));
t2s  = seconds(t2  - t1(1));
t3s  = seconds(t3  - t1(1));
t4s  = seconds(t4  - t1(1));
tPTs = seconds(tPT - t1(1));

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

tGrid  = t1s(keepIdx);
lc1Val = LC1.LC1(keepIdx);
lc2Val = lc2Val(keepIdx);
lc3Val = lc3Val(keepIdx);
lc4Val = lc4Val(keepIdx);
ptVal  = ptVal(keepIdx);

cellData = [lc1Val, lc2Val, lc3Val, lc4Val];
PT_vec   = ptVal .* (((2.52/2)^2)*pi);

fprintf('Matched rows: %d\n', sum(keepIdx));

mask     = ~any(cellData == 0, 2) & (PT_vec ~= 0);
cellData = cellData(mask, :);
PT_vec   = PT_vec(mask);
tGrid    = tGrid(mask);
fprintf('Rows with all nonzero readings: %d\n', sum(mask));

calMask  = tGrid >= 0 & tGrid <= 2000;
cellData = cellData(calMask, :);
PT_vec   = PT_vec(calMask);
tGrid    = tGrid(calMask);
fprintf('Rows in calibration window: %d\n', sum(calMask));

%% Save unfiltered calibration data
cellData_full = cellData;
PT_vec_full   = PT_vec;
tGrid_full    = tGrid;

%% Hysteresis check plot
dPT_check  = [0; diff(PT_vec_full)];
risingMask = dPT_check >= 0;

figure(20)
tiledlayout(2,2)
lcNames = {'LC1','LC2','LC3','LC4'};
for lc = 1:4
    nexttile
    scatter(PT_vec_full(risingMask),  cellData_full(risingMask,  lc), 3, 'b', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Rising')
    hold on
    scatter(PT_vec_full(~risingMask), cellData_full(~risingMask, lc), 3, 'r', 'filled', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Falling')
    xlabel('PT (lbf)'); ylabel('LC reading (lbf)')
    title(lcNames{lc}); legend
end
sgtitle('Hysteresis Check — Rising vs Falling Thrust')

%% Split by firing and fit per firing
dt_full   = diff(tGrid_full);
gapIdx    = find(dt_full > 1.0);
segStarts = [1; gapIdx+1];
segEnds   = [gapIdx; length(tGrid_full)];

fprintf('\nFound %d firings\n', length(segStarts));

firingCoeffs = zeros(length(segStarts), 5);
firingValid  = false(length(segStarts), 1);

for f = 1:length(segStarts)
    idx  = segStarts(f):segEnds(f);
    cd   = cellData_full(idx, :);
    pv   = PT_vec_full(idx);
    tSeg = tGrid_full(idx);

    % Rising only
    dPT_seg = [0; diff(pv)];
    rm      = dPT_seg >= 0;
    cd      = cd(rm, :);
    pv      = pv(rm);

    if length(pv) < 5
        fprintf('Firing %d: not enough data\n', f);
        continue
    end

    cb  = [cd, ones(size(cd,1),1)] \ pv;
    res = pv - (cd * cb(1:4) + cb(5));
    om  = abs(res - mean(res)) <= 3*std(res);
    cd  = cd(om,:); pv = pv(om);
    cb  = [cd, ones(size(cd,1),1)] \ pv;

    firingCoeffs(f,:) = cb';
    firingValid(f)    = true;

    est = cd * cb(1:4) + cb(5);
    fprintf('\n--- Firing %d (t=%.1f to %.1f s, %d rows) ---\n', ...
        f, tSeg(1), tSeg(end), length(pv));
    fprintf('LC1: %.4f  LC2: %.4f  LC3: %.4f  LC4: %.4f  Bias: %.2f\n', cb);
    fprintf('RMSE: %.2f  Max: %.2f  Mean: %.2f\n', ...
        sqrt(mean((pv-est).^2)), max(abs(pv-est)), mean(abs(pv-est)));
end

%% Plot LC pair sums vs PT
figure(5)
plot(tGrid_full, cellData_full(:,1) + cellData_full(:,2), 'b', 'DisplayName', 'LC1+LC2')
hold on
plot(tGrid_full, cellData_full(:,3) + cellData_full(:,4), 'r', 'DisplayName', 'LC3+LC4')
plot(tGrid_full, PT_vec_full, 'k--', 'DisplayName', 'PT')
legend; xlabel('Time (s)'); ylabel('Force (lbf)')
title('LC Pair Sums vs PT'); grid on

%% Test each firing's coefficients against full dataset
for f = 1:length(segStarts)
    if ~firingValid(f); continue; end
    fprintf('\n--- Testing Firing %d Coefficients ---\n', f);
    testCoeffs(firingCoeffs(f,:), 'integrated_1ms.csv')
end

%% Firing 4 self-test restricted to 0-7000 lbf operational range
idx  = segStarts(4):segEnds(4);
cd   = cellData_full(idx, :);
pv   = PT_vec_full(idx);
tSeg = tGrid_full(idx);

dPT_seg = [0; diff(pv) ./ diff(tSeg)];
rm      = dPT_seg >= 0 & pv >= 0 & pv <= 7000 & abs(dPT_seg) < 50;
cd      = cd(rm, :); pv = pv(rm); tSeg = tSeg(rm);

cb  = firingCoeffs(4,:);
est = cd * cb(1:4)' + cb(5);
err = est - pv;

fprintf('Firing 4 operational range (0-7000 lbf, rising, d1<50):\n');
fprintf('RMSE: %.2f lbf  Max: %.2f lbf  Mean: %.2f lbf  rows=%d\n', ...
    sqrt(mean(err.^2)), max(abs(err)), mean(abs(err)), length(pv));
fprintf('Mean error as %% of nominal (7000 lbf): %.2f%%\n', 100*mean(abs(err))/7000);
fprintf('RMSE as %% of nominal (7000 lbf): %.2f%%\n', 100*sqrt(mean(err.^2))/7000);

%% Find and plot the max error spike in Firing 4
idx  = segStarts(4):segEnds(4);
cd   = cellData_full(idx, :);
pv   = PT_vec_full(idx);
tSeg = tGrid_full(idx);

dPT_seg = [0; diff(pv) ./ diff(tSeg)];
rm      = dPT_seg >= 0 & pv >= 5000 & abs(dPT_seg) < 50;
cd      = cd(rm, :); pv = pv(rm); tSeg = tSeg(rm);

cb  = firingCoeffs(4,:);
est = cd * cb(1:4)' + cb(5);
err = est - pv;

[~, maxIdx] = max(abs(err));
fprintf('Max error spike: t=%.2f s  PT=%.1f lbf  Est=%.1f lbf  Err=%.1f lbf\n', ...
    tSeg(maxIdx), pv(maxIdx), est(maxIdx), err(maxIdx));

figure(15)
plot(tSeg, pv, 'k', 'DisplayName', 'PT')
hold on
plot(tSeg, est, 'r--', 'DisplayName', 'LC estimate')
plot(tSeg(maxIdx), pv(maxIdx), 'go', 'MarkerSize', 10, 'DisplayName', 'Max error')
legend; xlabel('Time (s)'); ylabel('Force (lbf)')
title('Firing 4 Steady State — Max Error Location')

%% Save best coefficients (Firing 4) to text file
bestFiring = 4;
cb = firingCoeffs(bestFiring, :);

fid = fopen('best_coefficients.txt', 'w');
fprintf(fid, 'Best Coefficients — Firing %d\n', bestFiring);
fprintf(fid, '=====================================\n');
fprintf(fid, 'LC1:  %.6f\n', cb(1));
fprintf(fid, 'LC2:  %.6f\n', cb(2));
fprintf(fid, 'LC3:  %.6f\n', cb(3));
fprintf(fid, 'LC4:  %.6f\n', cb(4));
fprintf(fid, 'Bias: %.4f lbf\n', cb(5));
fprintf(fid, '\n--- Operational Range (0-7000 lbf, rising, d1<50) ---\n');
fprintf(fid, 'RMSE:       %.2f lbf\n', sqrt(mean(err.^2)));
fprintf(fid, 'Max Error:  %.2f lbf\n', max(abs(err)));
fprintf(fid, 'Mean Error: %.2f lbf\n', mean(abs(err)));
fprintf(fid, 'Mean Error (%% of 7000 lbf nominal): %.2f%%\n', 100*mean(abs(err))/7000);
fprintf(fid, 'RMSE (%% of 7000 lbf nominal): %.2f%%\n', 100*sqrt(mean(err.^2))/7000);
fclose(fid);
fprintf('Coefficients saved to best_coefficients.txt\n');

%% Save validation plot
figure(16)
plot(tSeg, pv,  'k',  'DisplayName', 'PT (actual)')
hold on
plot(tSeg, est, 'r--','DisplayName', 'LC estimate')
legend; xlabel('Time (s)'); ylabel('Thrust (lbf)')
title(sprintf('Firing %d — Operational Range Validation (0-7000 lbf)', bestFiring))
grid on
saveas(figure(16), 'best_coefficients_validation.png')
fprintf('Validation plot saved to best_coefficients_validation.png\n');