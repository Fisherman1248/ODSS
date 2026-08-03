clc;
clear;
close all;

%% ===== Parameters =====
q_list = [1.5 2.0];

% Must be consistent with TX.m
Nscale = 6;

%% ===== Output folder =====
outDir = fullfile(pwd, 'ODSS_AWGN_FINAL_PDF');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ===== Style settings =====
% One color for each q
color_q = lines(length(q_list));

% One marker for each q
marker_list = {'o', 's', '^', 'd', 'v', '>', '<', 'p', 'h'};

fontName       = 'Times New Roman';
axisFontSize   = 14;
labelFontSize  = 16;
legendFontSize = 10.8;

simulationLineWidth = 1.75;
theoryLineWidth     = 1.50;
markerSize          = 6.0;

%% ===== Create figure =====
fig = figure('Color', 'w');

fig.Units = 'inches';
fig.Position = [1 1 7.2 5.0];

ax = axes(fig);

hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

hasAnyValidFile = false;
SNR_ref = [];

numExpectedCurves = 2*length(q_list);

legendHandles = gobjects(numExpectedCurves, 1);
legendLabels = cell(numExpectedCurves, 1);

legCount = 0;

%% ===== Load and plot each q result =====
for qq = 1:length(q_list)

    q = q_list(qq);

    %% Folder and receiver result filename
    folderName = sprintf('ODSS_AWGN_q%g_Nscale%d', q, Nscale);

    rxTag = sprintf('RX_q%g_Nscale%d', q, Nscale);

    rxFile = fullfile(folderName, [rxTag '.mat']);

    fprintf('Loading: %s\n', rxFile);

    %% Check file
    if exist(rxFile, 'file') ~= 2
        warning('File not found: %s', rxFile);
        continue;
    end

    %% Load receiver results
    S = load( ...
        rxFile, ...
        'BER_LS', ...
        'BER_theory', ...
        'SNR_dB');

    if ~isfield(S, 'BER_LS') || ~isfield(S, 'SNR_dB')
        warning('BER_LS or SNR_dB is missing from: %s', rxFile);
        continue;
    end

    BER_LS = S.BER_LS(:);
    SNR_dB = S.SNR_dB(:);

    if isfield(S, 'BER_theory')
        BER_theory = S.BER_theory(:);
    else
        BER_theory = [];
        warning('BER_theory is missing from: %s', rxFile);
    end

    %% Validate vector lengths
    if numel(BER_LS) ~= numel(SNR_dB)
        warning('BER_LS and SNR_dB have different lengths in: %s', ...
            rxFile);
        continue;
    end

    if ~isempty(BER_theory) && ...
            numel(BER_theory) ~= numel(SNR_dB)

        warning(['BER_theory and SNR_dB have different lengths ', ...
                 'in: %s'], rxFile);

        BER_theory = [];
    end

    %% Save reference SNR axis
    if isempty(SNR_ref)
        SNR_ref = SNR_dB;
    end

    %% Prevent zero values from disappearing on a logarithmic axis
    BER_LS_plot = BER_LS;
    BER_LS_plot(BER_LS_plot <= 0) = 1e-10;

    if ~isempty(BER_theory)
        BER_theory_plot = BER_theory;
        BER_theory_plot(BER_theory_plot <= 0) = 1e-10;
    end

    %% Plot settings
    curveColor = color_q(qq,:);

    markerIndex = mod(qq-1, length(marker_list)) + 1;
    curveMarker = marker_list{markerIndex};

    %% Simulation BER
    hSimulation = semilogy( ...
        ax, ...
        SNR_dB, ...
        BER_LS_plot, ...
        'Color', curveColor, ...
        'LineStyle', '-', ...
        'Marker', curveMarker, ...
        'LineWidth', simulationLineWidth, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', 'none', ...
        'MarkerIndices', 1:2:length(SNR_dB));

    legCount = legCount + 1;

    legendHandles(legCount) = hSimulation;
    legendLabels{legCount} = sprintf( ...
        '$q=%g$, simulation', q);

    %% Theoretical BER
    if ~isempty(BER_theory)

        hTheory = semilogy( ...
            ax, ...
            SNR_dB, ...
            BER_theory_plot, ...
            'Color', curveColor, ...
            'LineStyle', '--', ...
            'LineWidth', theoryLineWidth);

        legCount = legCount + 1;

        legendHandles(legCount) = hTheory;
        legendLabels{legCount} = sprintf( ...
            '$q=%g$, theory', q);
    end

    hasAnyValidFile = true;
end

%% ===== Stop if no valid result exists =====
if ~hasAnyValidFile
    close(fig);

    error(['No valid ODSS receiver files were found. ', ...
           'Check folderName and rxTag naming rules.']);
end

%% ===== Axes style =====
ax.FontName = fontName;
ax.FontSize = axisFontSize;
ax.LineWidth = 1.0;

ax.TickDir = 'in';
ax.TickLength = [0.012 0.012];

ax.YScale = 'log';

ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

ax.XGrid = 'on';
ax.YGrid = 'on';

ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

ax.GridAlpha = 0.25;

xlabel( ...
    ax, ...
    'Time-domain waveform SNR (dB)', ...
    'FontName', fontName, ...
    'FontSize', labelFontSize);

ylabel( ...
    ax, ...
    'Bit Error Rate', ...
    'FontName', fontName, ...
    'FontSize', labelFontSize);

% Leave the title empty for paper figures.
title(ax, '');

%% ===== Axis limits =====
ylim(ax, [1e-6 1]);

yticks(ax, [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1]);

yticklabels(ax, ...
    {'10^{-6}', '10^{-5}', '10^{-4}', ...
     '10^{-3}', '10^{-2}', '10^{-1}', '10^{0}'});

if ~isempty(SNR_ref)

    if numel(SNR_ref) >= 2 && max(SNR_ref) > min(SNR_ref)

        xlim(ax, [min(SNR_ref), max(SNR_ref)]);

    else
        xlim(ax, [SNR_ref(1)-1, SNR_ref(1)+1]);
    end
end

%% ===== Legend =====
legendHandles = legendHandles(1:legCount);
legendLabels = legendLabels(1:legCount);

lgd = legend( ...
    ax, ...
    legendHandles, ...
    legendLabels, ...
    'Interpreter', 'latex', ...
    'FontSize', legendFontSize, ...
    'Box', 'on', ...
    'Location', 'southwest');

lgd.NumColumns = 1;
lgd.ItemTokenSize = [24 9];

lgd.Color = 'white';
lgd.EdgeColor = [0.2 0.2 0.2];
lgd.LineWidth = 0.7;

%% ===== Export vector PDF =====
set(fig, 'InvertHardcopy', 'off');
set(fig, 'Color', 'w');

baseName = sprintf( ...
    'ODSS_AWGN_Nscale%d_BER', ...
    Nscale);

pdfFile = fullfile( ...
    outDir, ...
    [baseName '.pdf']);

exportgraphics( ...
    fig, ...
    pdfFile, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

fprintf('\nSaved PDF figure:\n%s\n', pdfFile);

%% ===== Also save MATLAB figure =====
figFile = fullfile( ...
    outDir, ...
    [baseName '.fig']);

savefig(fig, figFile);

fprintf('Saved MATLAB figure:\n%s\n', figFile);
fprintf('\nAll ODSS AWGN curves generated.\n');