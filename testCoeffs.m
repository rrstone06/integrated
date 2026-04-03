function testCoeffs(coeffs, csvFile)
% testCoeffs - Tests LC calibration coefficients against real PT data
% Usage: testCoeffs([c1, c2, c3, c4], 'integrated_1ms.csv')

%% Load and Format Data
data = readtable(csvFile);
data = data(:, [6:7, 13]);
data.Properties.VariableNames = {'time', 'value', 'name'};

LC1 = data(strcmp(data.name,'LC_1'), {'time','value'}); LC1.Properties.VariableNames{2} = 'LC1';
LC2 = data(strcmp(data.name,'LC_2'), {'time','value'}); LC2.Properties.VariableNames{2} = 'LC2';
LC3 = data(strcmp(data.name,'LC_3'), {'time','value'}); LC3.Properties.VariableNames{2} = 'LC3';
LC4 = data(strcmp(data.name,'LC_4'), {'time','value'}); LC4.Properties.VariableNames{2} = 'LC4';
PT  = data(strcmp(data.name,'PT_P_E_F_7'), {'time','value'}); PT.Properties.VariableNames{2} = 'PT';

%% Parse Timestamps
parseDT = @(c) datetime(regexprep(string(c), '(\d{2}:\d{2}:\d{2})Z', '$1.000Z'), ...
               'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss.SSSZ', 'TimeZone', 'UTC');

t1  = seconds(parseDT(LC1.time) - parseDT(LC1.time(1)));
t2  = seconds(parseDT(LC2.time) - parseDT(LC1.time(1)));
t3  = seconds(parseDT(LC3.time) - parseDT(LC1.time(1)));
t4  = seconds(parseDT(LC4.time) - parseDT(LC1.time(1)));
tPT = seconds(parseDT(PT.time)  - parseDT(LC1.time(1)));

%% Interpolate to 1ms Grid
tStart = max([t1(1), t2(1), t3(1), t4(1), tPT(1)]);
tEnd   = min([t1(end), t2(end), t3(end), t4(end), tPT(end)]);
tGrid  = (tStart:0.001:tEnd)';

lc1i = interp1(t1,  LC1.LC1, tGrid, 'linear');
lc2i = interp1(t2,  LC2.LC2, tGrid, 'linear');
lc3i = interp1(t3,  LC3.LC3, tGrid, 'linear');
lc4i = interp1(t4,  LC4.LC4, tGrid, 'linear');
pti  = interp1(tPT, PT.PT,   tGrid, 'linear') .* 4.98;

%% Remove Zero Rows
lcMat = [lc1i, lc2i, lc3i, lc4i];
mask  = ~any(lcMat == 0, 2) & (pti ~= 0);
lcMat = lcMat(mask, :);
pti   = pti(mask);
tGrid = tGrid(mask);

%% Apply Coefficients
forceEst = lcMat * coeffs(:);

%% Compute Error
err    = forceEst - pti;
rmse   = sqrt(mean(err.^2));
maxErr = max(abs(err));

fprintf('RMSE:      %.4f\n', rmse);
fprintf('Max Error: %.4f\n', maxErr);

figure(3)
plot(tGrid, pti,      'k',  'DisplayName', 'PT (actual, lbf)')
hold on
plot(tGrid, forceEst, 'r--','DisplayName', 'LC estimate (lbf)')
legend
xlabel('Time (s)')
ylabel('Thrust (lbf)')
title('Coefficient Validation')

figure(4)
plot(tGrid, err)
xlabel('Time (s)')
ylabel('Error (lbf)')
title('Residual Error')
yline(0, 'k--')

end