%% Reproduce Fig. 2(b)-style band dispersion from COMSOL eigenfrequency scan
% Data columns:
% 1 delta, 2 k ratio, 3 complex eigenfrequency, 4 Q, 5-8 Stokes S0-S3.
% The paired degenerate modes are taken as the 5th and 6th eigen-solutions
% within each fixed (delta, k) scan point.

clear; clc; close all;

xlsxPath = "E:\Pang\Project\2025_11\Degeneracy\Data\S2n Symmetry\S6 Symmetry\Band Structure\Continuous\Nonreciprocity,delta=0&0.08,k=(-0.06,0.06,length=0.005).xlsx";
outDir = fullfile(pwd, "outputs");
if ~exist(outDir, "dir")
    mkdir(outDir);
end

modePair = [5, 6];       % 1-based eigenmode indices in each (delta, k) block
delta0 = 0;
deltaB = 0.08;

raw = readcell(xlsxPath);
raw = raw(~all(cellfun(@isEmptyCell, raw), 2), :);

delta = cellfun(@toDouble, raw(:, 1));
k = cellfun(@toDouble, raw(:, 2));
freqComplex = cellfun(@parseComplex, raw(:, 3));
qFactor = cellfun(@toDouble, raw(:, 4));
s0 = cellfun(@toDouble, raw(:, 5));
s1 = cellfun(@toDouble, raw(:, 6));
s2 = cellfun(@toDouble, raw(:, 7));
s3 = cellfun(@toDouble, raw(:, 8));

valid = ~isnan(delta) & ~isnan(k) & ~isnan(real(freqComplex));
T = table(delta(valid), k(valid), freqComplex(valid), qFactor(valid), ...
    s0(valid), s1(valid), s2(valid), s3(valid), find(valid), ...
    'VariableNames', {'delta','k','freq','Q','S0','S1','S2','S3','sourceRow'});
T.freqReal = real(T.freq);

T = sortrows(T, {'delta','k','sourceRow'});
[G, ~] = findgroups(T.delta, T.k);
T.modeIndex = zeros(height(T), 1);
for g = 1:max(G)
    rows = find(G == g);
    T.modeIndex(rows) = (1:numel(rows)).';
end

band0 = extractBandPair(T, delta0, modePair);
bandB = extractBandPair(T, deltaB, modePair);

fprintf("delta = %.3g, k = 0: mode %d = %.6f, mode %d = %.6f\n", ...
    delta0, modePair(1), gammaValue(band0.k, band0.lower), ...
    modePair(2), gammaValue(band0.k, band0.upper));
fprintf("delta = %.3g, k = 0: lower = %.6f, upper = %.6f, splitting = %.6f\n", ...
    deltaB, gammaValue(bandB.k, bandB.lower), gammaValue(bandB.k, bandB.upper), ...
    gammaValue(bandB.k, bandB.upper) - gammaValue(bandB.k, bandB.lower));

fig = figure('Color', 'w', 'Units', 'centimeters', 'Position', [4 4 8.5 8.0]);
ax = axes(fig);
ax.Position = [0.18 0.42 0.76 0.50];
hold(ax, 'on');
try
    ax.Toolbar.Visible = 'off';
catch
end
try
    disableDefaultInteractivity(ax);
catch
end

% Paper-style colors: black for delta = 0, red/blue for the lifted pair.
black = [0.10, 0.10, 0.10];
red = [0.95, 0.05, 0.05];
blue = [0.05, 0.15, 0.95];
gray = [0.45, 0.45, 0.45];

hDelta0 = plot(ax, band0.k, band0.upper, '-', 'Color', black, 'LineWidth', 1.05);
plot(ax, band0.k, band0.lower, '-', 'Color', black, 'LineWidth', 1.05, ...
    'HandleVisibility', 'off');
hUpper = plot(ax, bandB.k, bandB.upper, '-', 'Color', red, 'LineWidth', 1.35);
hLower = plot(ax, bandB.k, bandB.lower, '-', 'Color', blue, 'LineWidth', 1.35);

kGamma = 0;
f0Gamma = mean([gammaValue(band0.k, band0.lower), gammaValue(band0.k, band0.upper)]);
fUpperGamma = gammaValue(bandB.k, bandB.upper);
fLowerGamma = gammaValue(bandB.k, bandB.lower);

plot(ax, kGamma, f0Gamma, 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', gray, 'MarkerEdgeColor', gray);
plot(ax, kGamma, fUpperGamma, 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', red, 'MarkerEdgeColor', red);
plot(ax, kGamma, fLowerGamma, 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', blue, 'MarkerEdgeColor', blue);

xline(ax, kGamma, '--', 'Color', [0.68, 0.68, 0.68], 'LineWidth', 0.9);

xmin = min([band0.k; bandB.k]);
xmax = max([band0.k; bandB.k]);
xlim(ax, [xmin, xmax]);
ylimPad = 0.08 * range([band0.lower; band0.upper; bandB.lower; bandB.upper]);
ymin = min([band0.lower; band0.upper; bandB.lower; bandB.upper]) - ylimPad;
ymax = max([band0.lower; band0.upper; bandB.lower; bandB.upper]) + ylimPad;
ylim(ax, [ymin, ymax]);

box(ax, 'on');
ax.LineWidth = 0.9;
ax.FontName = 'Arial';
ax.FontSize = 8;
ax.TickDir = 'in';
ax.Layer = 'top';
ax.XTick = [xmin, 0, xmax];
ax.XTickLabel = {'M', '\Gamma', 'K'};

xlabel(ax, 'k (2\pi/a)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'Re(\omega) (THz)', 'FontName', 'Arial', 'FontSize', 9);

leg = legend(ax, [hDelta0, hUpper, hLower], ...
    {'\delta = 0', 'Upper band, \delta = 0.08', 'Lower band, \delta = 0.08'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', 2);
leg.Units = 'normalized';
leg.Position = [0.12 0.11 0.82 0.14];
leg.FontSize = 7;

text(ax, 0.02, 0.94, '(b)', 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);

pngPath = fullfile(outDir, "fig2b_upper_lower_bands_matlab.png");
pdfPath = fullfile(outDir, "fig2b_upper_lower_bands_matlab.pdf");
exportgraphics(fig, pngPath, 'Resolution', 600);
set(fig, 'Renderer', 'painters');
print(fig, pdfPath, '-dpdf', '-painters');

fprintf("Saved PNG: %s\n", pngPath);
fprintf("Saved PDF: %s\n", pdfPath);

function band = extractBandPair(T, deltaValue, modePair)
    tol = 1e-9;
    sub = T(abs(T.delta - deltaValue) < tol & ismember(T.modeIndex, modePair), :);
    if isempty(sub)
        error("No rows found for delta = %.6g and modes [%d, %d].", ...
            deltaValue, modePair(1), modePair(2));
    end

    kVals = unique(sub.k, 'sorted');
    lower = nan(size(kVals));
    upper = nan(size(kVals));
    for ii = 1:numel(kVals)
        rows = sub(abs(sub.k - kVals(ii)) < tol, :);
        if height(rows) ~= 2
            error("Expected two selected modes at delta = %.6g, k = %.6g, found %d.", ...
                deltaValue, kVals(ii), height(rows));
        end
        pair = sort(rows.freqReal);
        lower(ii) = pair(1);
        upper(ii) = pair(2);
    end

    band = table(kVals, lower, upper, ...
        'VariableNames', {'k','lower','upper'});
end

function y0 = gammaValue(k, y)
    [dmin, idx] = min(abs(k));
    if dmin > 1e-10
        error("No Gamma-point k = 0 row was found.");
    end
    y0 = y(idx);
end

function value = toDouble(x)
    if isnumeric(x)
        value = double(x);
    elseif ismissing(string(x)) || strlength(strtrim(string(x))) == 0
        value = NaN;
    else
        value = str2double(string(x));
    end
end

function z = parseComplex(x)
    if isnumeric(x)
        z = complex(double(x), 0);
        return;
    end

    s = erase(strtrim(string(x)), " ");
    s = replace(s, "I", "i");
    s = replace(s, "j", "i");

    if strlength(s) == 0 || ismissing(s)
        z = complex(NaN, NaN);
        return;
    end

    z0 = str2double(s);
    if ~isnan(real(z0))
        z = z0;
        return;
    end

    token = regexp(s, '^([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)([+-]\d*\.?\d+(?:[eE][+-]?\d+)?)i$', ...
        'tokens', 'once');
    if isempty(token)
        z = complex(NaN, NaN);
    else
        z = complex(str2double(token{1}), str2double(token{2}));
    end
end

function tf = isEmptyCell(x)
    tf = isempty(x) || (isstring(x) && ismissing(x)) || ...
        (ischar(x) && isempty(strtrim(x))) || ...
        (isnumeric(x) && isscalar(x) && isnan(x));
end
