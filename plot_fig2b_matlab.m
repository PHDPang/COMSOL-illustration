%% 根据 COMSOL 本征频率扫描结果绘制复本征频率的实部和虚部
% 数据列：1-delta，2-k/ratio，3-复本征频率，4-Q，5~8-S0~S3。
% 对每个固定的 (delta, k) 扫描点，选取第 5、6 个本征解作为目标模式。

clear; clc; close all;

xlsxPath = "E:\Pang\Project\2025_11\Degeneracy\Data\S2n Symmetry\S6 Symmetry\Band Structure\Continuous\Nonreciprocity,delta=0&0.08,k=(-0.04,0.04,length=0.005).xlsx";
outDir = fullfile(pwd, "outputs");
% 若输出目录不存在则自动创建。
if ~exist(outDir, "dir")
    mkdir(outDir);
end

modePair = [5, 6];       % 每个 (delta, k) 数据块内的本征模式序号
delta0 = 0;
deltaB = 0.08;
c_const = 299792458;     % 真空光速，单位：m/s
lattice = 840e-9;        % 晶格常数，单位：m
f0 = c_const / lattice / 1e12;  % 归一化频率，单位：THz

% 读取原始单元格，并删除整行为空的数据。
raw = readcell(xlsxPath);
raw = raw(~all(cellfun(@isEmptyCell, raw), 2), :);

% 将各列转换为数值；复本征频率单独解析为复数。
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
% 消除 Excel 浮点存储误差，确保理论上相同的 k 被归入同一组。
T.k = round(T.k, 12);

% 保留源文件中的本征解顺序，并在每个 (delta, k) 组内重新编号。
T = sortrows(T, {'delta','k','sourceRow'});
[G, ~] = findgroups(T.delta, T.k);
T.modeIndex = zeros(height(T), 1);
for g = 1:max(G)
    rows = find(G == g);
    T.modeIndex(rows) = (1:numel(rows)).';
end

% 分别提取 delta=0 和 delta=0.08 时的上下能带。
band0 = extractBandPair(T, delta0, modePair);
bandB = extractBandPair(T, deltaB, modePair);

fprintf("delta = %.3g, k = 0: mode %d = %.6f, mode %d = %.6f\n", ...
    delta0, modePair(1), gammaValue(band0.k, band0.lower), ...
    modePair(2), gammaValue(band0.k, band0.upper));
fprintf("delta = %.3g, k = 0: lower = %.6f, upper = %.6f, splitting = %.6f\n", ...
    deltaB, gammaValue(bandB.k, bandB.lower), gammaValue(bandB.k, bandB.upper), ...
    gammaValue(bandB.k, bandB.upper) - gammaValue(bandB.k, bandB.lower));

% 保留本征频率虚部的原始 THz 数据，不进行归一化。
imagBand0 = table(band0.k, band0.lowerImag, band0.upperImag, ...
    'VariableNames', {'k','lower','upper'});
imagBandB = table(bandB.k, bandB.lowerImag, bandB.upperImag, ...
    'VariableNames', {'k','lower','upper'});

% 将频率除以 f0，后续绘图统一使用无量纲频率。
band0.lower = band0.lower / f0;
band0.upper = band0.upper / f0;
bandB.lower = bandB.lower / f0;
bandB.upper = bandB.upper / f0;

% 仅对绘图曲线进行保形插值，避免改变原始数据及其物理含义。
interpFactor = 10;
plotBand0 = interpolateBandForPlot(band0, interpFactor);
plotBandB = interpolateBandForPlot(bandB, interpFactor);
plotImagBand0 = interpolateBandForPlot(imagBand0, interpFactor);
plotImagBandB = interpolateBandForPlot(imagBandB, interpFactor);

fig = figure(1);
set(fig, 'Name', 'Real part (a)', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [4 4 8.5 8.0]);
ax = axes(fig);
ax.Position = [0.18 0.42 0.76 0.50];
hold(ax, 'on');
% 隐藏交互工具栏，避免导出图片时出现额外控件。
try
    ax.Toolbar.Visible = 'off';
catch
end
try
    disableDefaultInteractivity(ax);
catch
end

% 参照论文配色：delta=0 为黑色，解除简并后的上下支分别为红、蓝色。
black = [0.10, 0.10, 0.10];
red = [0.95, 0.05, 0.05];
blue = [0.05, 0.15, 0.95];
gray = [0.45, 0.45, 0.45];

hDelta0 = plot(ax, plotBand0.k, plotBand0.upper, '-', 'Color', black, 'LineWidth', 1.05);
plot(ax, plotBand0.k, plotBand0.lower, '-', 'Color', black, 'LineWidth', 1.05, ...
    'HandleVisibility', 'off');
hUpper = plot(ax, plotBandB.k, plotBandB.upper, '-', 'Color', red, 'LineWidth', 1.35);
hLower = plot(ax, plotBandB.k, plotBandB.lower, '-', 'Color', blue, 'LineWidth', 1.35);

% 提取并标记 Gamma 点处的简并位置及分裂后的上下支。
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

% 根据全部能带数据统一设置显示范围，并预留少量纵向边距。
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

% 横坐标每隔 0.02 设置刻度，并将高对称点标记放在数值下一行。
xTicks = xmin:0.02:xmax;
xTickLabels = arrayfun(@(x) sprintf('%.2f', x), xTicks, ...
    'UniformOutput', false);
[~, idxM] = min(abs(xTicks - xmin));
[~, idxGamma] = min(abs(xTicks));
[~, idxK] = min(abs(xTicks - xmax));
xTickLabels{idxM} = sprintf('%.2f\\newlineM', xTicks(idxM));
xTickLabels{idxGamma} = sprintf('%.2f\\newline\\Gamma', xTicks(idxGamma));
xTickLabels{idxK} = sprintf('%.2f\\newlineK', xTicks(idxK));
ax.XTick = xTicks;
ax.XTickLabel = xTickLabels;
ax.TickLabelInterpreter = 'tex';

xlabel(ax, 'k (2\pi/a)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax, 'f/f_0', 'FontName', 'Arial', 'FontSize', 9);

leg = legend(ax, [hDelta0, hUpper, hLower], ...
    {'\delta = 0', 'Upper band, \delta = 0.08', 'Lower band, \delta = 0.08'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', 2);
leg.Units = 'normalized';
leg.Position = [0.12 0.11 0.82 0.14];
leg.FontSize = 7;

text(ax, 0.02, 0.94, '(a)', 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);

% 同时导出高分辨率 PNG 和矢量 PDF。
pngPath = fullfile(outDir, "fig2a_real_upper_lower_bands_matlab.png");
pdfPath = fullfile(outDir, "fig2a_real_upper_lower_bands_matlab.pdf");
exportgraphics(fig, pngPath, 'Resolution', 600);
set(fig, 'Renderer', 'painters');
print(fig, pdfPath, '-dpdf', '-painters');

fprintf("Saved PNG: %s\n", pngPath);
fprintf("Saved PDF: %s\n", pdfPath);

% 绘制复本征频率虚部，样式与实部图保持一致。
figImag = figure(2);
set(figImag, 'Name', 'Imaginary part (b)', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [13 4 8.5 8.0]);
axImag = axes(figImag);
axImag.Position = [0.18 0.42 0.76 0.50];
hold(axImag, 'on');
try
    axImag.Toolbar.Visible = 'off';
catch
end
try
    disableDefaultInteractivity(axImag);
catch
end

hImagDelta0 = plot(axImag, plotImagBand0.k, plotImagBand0.upper, '-', ...
    'Color', black, 'LineWidth', 1.05);
plot(axImag, plotImagBand0.k, plotImagBand0.lower, '-', ...
    'Color', black, 'LineWidth', 1.05, 'HandleVisibility', 'off');
hImagUpper = plot(axImag, plotImagBandB.k, plotImagBandB.upper, '-', ...
    'Color', red, 'LineWidth', 1.35);
hImagLower = plot(axImag, plotImagBandB.k, plotImagBandB.lower, '-', ...
    'Color', blue, 'LineWidth', 1.35);

% 标记 Gamma 点处各模式的原始虚部。
plot(axImag, kGamma, gammaValue(imagBand0.k, imagBand0.upper), 'o', ...
    'MarkerSize', 5.5, 'MarkerFaceColor', gray, 'MarkerEdgeColor', gray);
plot(axImag, kGamma, gammaValue(imagBand0.k, imagBand0.lower), 'o', ...
    'MarkerSize', 5.5, 'MarkerFaceColor', gray, 'MarkerEdgeColor', gray);
plot(axImag, kGamma, gammaValue(imagBandB.k, imagBandB.upper), 'o', ...
    'MarkerSize', 5.5, 'MarkerFaceColor', red, 'MarkerEdgeColor', red);
plot(axImag, kGamma, gammaValue(imagBandB.k, imagBandB.lower), 'o', ...
    'MarkerSize', 5.5, 'MarkerFaceColor', blue, 'MarkerEdgeColor', blue);
xline(axImag, kGamma, '--', 'Color', [0.68, 0.68, 0.68], 'LineWidth', 0.9);

xlim(axImag, [xmin, xmax]);
imagValues = [imagBand0.lower; imagBand0.upper; imagBandB.lower; imagBandB.upper];
imagPad = 0.08 * range(imagValues);
if imagPad == 0
    imagPad = max(abs(imagValues(1)) * 0.08, eps);
end
ylim(axImag, [min(imagValues) - imagPad, max(imagValues) + imagPad]);

box(axImag, 'on');
axImag.LineWidth = 0.9;
axImag.FontName = 'Arial';
axImag.FontSize = 8;
axImag.TickDir = 'in';
axImag.Layer = 'top';
axImag.XTick = xTicks;
axImag.XTickLabel = xTickLabels;
axImag.TickLabelInterpreter = 'tex';

xlabel(axImag, 'k (2\pi/a)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(axImag, 'Im(f) (THz)', 'FontName', 'Arial', 'FontSize', 9);

legImag = legend(axImag, [hImagDelta0, hImagUpper, hImagLower], ...
    {'\delta = 0', 'Upper band, \delta = 0.08', 'Lower band, \delta = 0.08'}, ...
    'Orientation', 'horizontal', 'Box', 'off', 'NumColumns', 2);
legImag.Units = 'normalized';
legImag.Position = [0.12 0.11 0.82 0.14];
legImag.FontSize = 7;

text(axImag, 0.02, 0.94, '(b)', 'Units', 'normalized', ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);

pngImagPath = fullfile(outDir, "fig2b_imag_upper_lower_bands_matlab.png");
pdfImagPath = fullfile(outDir, "fig2b_imag_upper_lower_bands_matlab.pdf");
exportgraphics(figImag, pngImagPath, 'Resolution', 600);
set(figImag, 'Renderer', 'painters');
print(figImag, pdfImagPath, '-dpdf', '-painters');

fprintf("Saved PNG: %s\n", pngImagPath);
fprintf("Saved PDF: %s\n", pdfImagPath);

function plotBand = interpolateBandForPlot(band, interpFactor)
    % 使用保形分段三次插值加密横坐标，减少折线感并抑制非物理过冲。
    pointCount = (height(band) - 1) * interpFactor + 1;
    kFine = linspace(band.k(1), band.k(end), pointCount).';
    lowerFine = interp1(band.k, band.lower, kFine, 'pchip');
    upperFine = interp1(band.k, band.upper, kFine, 'pchip');
    plotBand = table(kFine, lowerFine, upperFine, ...
        'VariableNames', {'k','lower','upper'});
end

function band = extractBandPair(T, deltaValue, modePair)
    % 按实部划分上下支，并同步保留每条能带对应的虚部。
    tol = 1e-9;
    sub = T(abs(T.delta - deltaValue) < tol & ismember(T.modeIndex, modePair), :);
    if isempty(sub)
        error("No rows found for delta = %.6g and modes [%d, %d].", ...
            deltaValue, modePair(1), modePair(2));
    end

    kVals = unique(sub.k, 'sorted');
    lower = nan(size(kVals));
    upper = nan(size(kVals));
    lowerImag = nan(size(kVals));
    upperImag = nan(size(kVals));
    for ii = 1:numel(kVals)
        rows = sub(abs(sub.k - kVals(ii)) < tol, :);
        if height(rows) ~= 2
            error("Expected two selected modes at delta = %.6g, k = %.6g, found %d.", ...
                deltaValue, kVals(ii), height(rows));
        end
        [pairReal, order] = sort(rows.freqReal);
        pairFreq = rows.freq(order);
        lower(ii) = pairReal(1);
        upper(ii) = pairReal(2);
        lowerImag(ii) = imag(pairFreq(1));
        upperImag(ii) = imag(pairFreq(2));
    end

    band = table(kVals, lower, upper, lowerImag, upperImag, ...
        'VariableNames', {'k','lower','upper','lowerImag','upperImag'});
end

function y0 = gammaValue(k, y)
    % 返回最接近 k=0 的 Gamma 点数值，并检查该点确实存在。
    [dmin, idx] = min(abs(k));
    if dmin > 1e-10
        error("No Gamma-point k = 0 row was found.");
    end
    y0 = y(idx);
end

function value = toDouble(x)
    % 兼容 Excel 中的数值单元格、文本数值和空单元格。
    if isnumeric(x)
        value = double(x);
    elseif ismissing(string(x)) || strlength(strtrim(string(x))) == 0
        value = NaN;
    else
        value = str2double(string(x));
    end
end

function z = parseComplex(x)
    % 将 COMSOL 导出的复数字符串统一转换为 MATLAB 复数。
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

    % str2double 失败时，用正则表达式解析“实部+虚部i”的形式。
    token = regexp(s, '^([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)([+-]\d*\.?\d+(?:[eE][+-]?\d+)?)i$', ...
        'tokens', 'once');
    if isempty(token)
        z = complex(NaN, NaN);
    else
        z = complex(str2double(token{1}), str2double(token{2}));
    end
end

function tf = isEmptyCell(x)
    % 判断不同数据类型下的空单元格。
    tf = isempty(x) || (isstring(x) && ismissing(x)) || ...
        (ischar(x) && isempty(strtrim(x))) || ...
        (isnumeric(x) && isscalar(x) && isnan(x));
end
