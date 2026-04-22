function testCoeffsPiecewise(binCoeffsAll, binValid, binEdges, csvFile)

%% Load and Format Data
data = readtable(csvFile);
data = data(:, [6:7, 13]);
data.Properties.VariableNames = {'time', 'value', 'name'};

LC1 = data(strcmp(data.name,'LC_1'), {'time','value'}); LC1.Properties.VariableNames{2} = 'LC1';
LC2 = data(strcmp(data.name,'LC_2'), {'time','value'}); LC2.Properties.VariableNames{2} = 'LC2';
LC3 = data(strcmp(data.name,'LC_3'), {'time','value'}); LC3.Properties.VariableNames{2} = 'LC3';
LC4 = data(strcmp(data.name,'LC_4'), {'time','value'}); LC4.Properties.VariableNames{2} = 'LC4';
PT  = data(strcmp(data.name,'PT_P_E_F_7'), {'time','value'}); PT.Properties.VariableNames{2} = 'PT';

parseDT = @(c) datetime(regexprep(string(c), '(\d{2}:\d{2}:\d{2})Z', '$1.000Z'), ...
    'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSSZ', 'TimeZone', 'UTC');

t1  = seconds(parseDT(LC1.time) - parseDT(LC1.time(1)));
t2  = seconds(parseDT(LC2.time) - parseDT(LC1.time(1)));
t3  = seconds(parseDT(LC3.time) - parseDT(LC1.time(1)));
t4  = seconds(parseDT(LC4.time) - parseDT(LC1.time(1)));
tPT = seconds(parseDT(PT.time)  - parseDT(LC1.time(1)));

tStart = max([t1(1), t2(1), t3(1), t4(1), tPT(1)]);
tEnd   = min([t1(end), t2(end), t3(end), t4(end), tPT(end), 2000]);
tGrid  = (tStart:0.001:tEnd)';

lc1i = interp1(t1,  LC1.LC1, tGrid, 'linear');
lc2i = interp1(t2,  LC2.LC2, tGrid, 'linear');
lc3i = interp1(t3,  LC3.LC3, tGrid, 'linear');
lc4i = interp1(t4,  LC4.LC4, tGrid, 'linear');
pti  = interp1(tPT, PT.PT,   tGrid, 'linear') .* (((2.52/2)^2)*pi);

lcMat = [lc1i, lc2i, lc3i, lc4i];
mask  = ~any(lcMat == 0, 2) & (pti ~= 0);
lcMat = lcMat(mask, :);
pti   = pti(mask);
tGrid = tGrid(mask);

% Rising only
dPT        = [0; diff(pti)];
risingMask = dPT >= 0;
lcMat      = lcMat(risingMask, :);
pti        = pti(risingMask);
tGrid      = tGrid(risingMask);

%% Apply piecewise coefficients
forceEst = zeros(size(pti));
nBins    = length(binEdges) - 1;

for b = 1:nBins
    if ~binValid(b); continue; end
    lo = binEdges(b);
    hi = binEdges(b+1);
    bm = pti >= lo & pti < hi;
    if sum(bm) == 0; continue; end
    cb = binCoeffsAll(b,:);
    forceEst(bm) = lcMat(bm,:) * cb(1:4)' + cb(5);
end

% Fill gaps with nearest valid bin
unfilled = forceEst == 0 & pti > 0;
if any(unfilled)
    validBins = find(binValid);
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;
    for i = find(unfilled)'
        [~, nearest] = min(abs(binCenters(validBins) - pti(i)));
        b  = validBins(nearest);
        cb = binCoeffsAll(b,:);
        forceEst(i) = lcMat(i,:) * cb(1:4)' + cb(5);
    end
end

%% Error metrics
err    = forceEst - pti;
rmse   = sqrt(mean(err.^2));
maxErr = max(abs(err));

fprintf('RMSE:       %.4f lbf\n', rmse);
fprintf('Max Error:  %.4f lbf\n', maxErr);
fprintf('Mean Error: %.4f lbf\n', mean(abs(err)));

threshPcts = [1, 2, 5, 10];
fprintf('\n--- Points within x%% of max error (%.1f lbf) ---\n', maxErr);
for t = threshPcts
    tol_lbf = t/100 * maxErr;
    pct_in  = 100 * mean(abs(err) <= tol_lbf);
    fprintf('  Within %d%% of max error (%.0f lbf): %.1f%% of points\n', t, tol_lbf, pct_in);
end

cullMask = abs(err) <= 0.95 * maxErr;
errCull  = err(cullMask);
fprintf('\n--- After culling top 5%% error points ---\n');
fprintf('RMSE:       %.4f lbf\n', sqrt(mean(errCull.^2)));
fprintf('Max Error:  %.4f lbf\n', max(abs(errCull)));
fprintf('Mean Error: %.4f lbf\n', mean(abs(errCull)));

%% Plots
figure(6)
plot(tGrid, pti,      'k',  'DisplayName', 'PT (actual)')
hold on
plot(tGrid, forceEst, 'r--','DisplayName', 'LC piecewise estimate')
legend; xlabel('Time (s)'); ylabel('Thrust (lbf)')
title('Piecewise Coefficient Validation')

figure(7)
plot(tGrid, err)
yline(0,'k--'); xlabel('Time (s)'); ylabel('Error (lbf)')
title('Piecewise Residual Error')

end