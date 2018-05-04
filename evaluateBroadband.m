function [out] = evaluateBroadband(spikeRate, broadband, params)

% Various calculations that reflects quality of recovery by broadband estimation
% relative to noiseless time course, e.g. regression

t = params.simulation.t/params.simulation.srate;

% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Extracted broadband, averaged across trials
%meanBroadband = geomean(broadband,2);
meanBroadband = mean(broadband,2);

out = [];

% Regression
stats = regstats(spikeRate(idx),meanBroadband(idx));
out.regress.rsq       = stats.rsquare;
out.regress.residuals = stats.r;
out.regress.sse       = stats.fstat.sse;

% SNR per band (?)
out.snr = [];

% Deviation from input
out.deviation = [];

% Linear/non-linear fits?
out.fit = [];

% Plot?
switch params.plot.on
    case 'yes'
        
        fH = figure;  %set(fH, 'Color', 'w');
        
        % Subtract 'prestim' baseline
        baseline = meanBroadband(t > -1 & t < 0);
        meanBroadband = meanBroadband(idx) - mean(baseline);
        %meanBroadband = meanBroadband(idx);
        
        % Scale for plotting
        mnToPlot = meanBroadband / norm(meanBroadband);
        spikeRateToPlot = spikeRate(idx) / norm(spikeRate(idx));
        
        plot(t(idx), mnToPlot, t(idx), spikeRateToPlot, 'k-', 'LineWidth', params.plot.lnwdth)
        title(params.analysis.methodstr)
        set(gca, 'XLim', params.plot.xl);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)
        xlabel('Time (s)')
        ylabel('Response')
end