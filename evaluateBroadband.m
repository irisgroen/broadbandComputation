function [results, params] = evaluateBroadband(spikeRate, broadband, params)

% Various calculations that reflects quality of recovery by broadband estimation
% relative to noiseless time course, e.g. regression

t = params.simulation.t/params.simulation.srate;

% Get 'calibrated' responses to scale estimatedBroadband with
[params] = calibrateResponseLevel(params);

% Extracted broadband, averaged across trials
%meanBroadband = geomean(broadband,2);
meanBroadband = mean(broadband,2);

% Scale broadband to calibrated response
meanBroadband = meanBroadband / (params.analysis.calibration(2)-params.analysis.calibration(1));

% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Compare idealized rate with the estimated broadband
results = [];

% Regression
stats = regstats(spikeRate(idx),meanBroadband(idx));
results.regress.rsq       = stats.rsquare;
results.regress.residuals = stats.r;
results.regress.sse       = stats.fstat.sse;

% SNR per band (?)
results.snr = [];

% Deviation from input
results.deviation = [];

% Linear/non-linear fits?
results.fit = [];

% Plot?
switch params.plot.on
    case 'yes'
        
        fH = figure;  set(fH, 'Color', 'w');
        
        % Subtract 'prestim' baseline
        baseline = meanBroadband(t > -1 & t < 0);
        meanBroadband = meanBroadband(idx) - mean(baseline);
        
        % Scale for plotting
        mnToPlot = meanBroadband; %/ norm(meanBroadband);
        spikeRateToPlot = spikeRate(idx); %/ norm(spikeRate(idx));
        
        plot(t(idx), spikeRateToPlot, t(idx), mnToPlot, 'k-', 'LineWidth', params.plot.lnwdth)
        title(params.analysis.methodstr)
        set(gca, 'XLim', params.plot.xl);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)
        xlabel('Time (s)')
        ylabel('Response')
        legend({'Noiseless time series', 'Estimated broadband signal'}, 'Location', 'NorthWest')

end