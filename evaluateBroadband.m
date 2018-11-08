function [results, params] = evaluateBroadband(spikeRate, estimatedBroadband, params)

% Various calculations that reflects quality of recovery by broadband estimation
% relative to noiseless time course, e.g. regression

t = params.simulation.t/params.simulation.srate;

% Average extracted broadband across trials
%meanBroadband = geomean(estimatedBroadband,2);
meanBroadband = mean(estimatedBroadband,2);

% Get 'calibrated' responses to scale estimatedBroadband with
[params] = calibrateResponseLevel(params);

% Scale broadband to calibrated response using calibrate function
meanBroadbandCalibrated = params.analysis.calibrate(meanBroadband);

% Clip time series to avoid edge artifacts
idx = t > -0.5 & t < 1;

% Compare idealized rate with the estimated broadband
results = [];

% Regression
stats = regstats(spikeRate(idx),meanBroadbandCalibrated(idx));
results.regress.rsq       = stats.rsquare;
results.regress.residuals = stats.r;
results.regress.sse       = stats.fstat.sse;

% Other measures? e.g SNR per band, deviation from input

% Plot
switch params.plot.on
    case 'yes'
        
        fH = figure;  set(fH, 'Color', 'w');
        
        % Subtract 'prestim' baseline
        %baseline = meanBroadbandCalibrated(t > -1 & t < 0);
        %meanBroadbandCalibrated = meanBroadbandCalibrated(idx) - mean(baseline);
        
        % Scale for plotting
        % mnToPlot = meanBroadband / norm(meanBroadband);
        % spikeRateToPlot = spikeRate(idx); %/ norm(spikeRate(idx));
        
        % Clip
        mnToPlot = meanBroadbandCalibrated(idx);
        spikeRateToPlot = spikeRate(idx);
        
        plot(t(idx), spikeRateToPlot, t(idx), mnToPlot, 'k-', 'LineWidth', params.plot.lnwdth)
        title(params.analysis.methodstr)
        set(gca, 'XLim', params.plot.xl);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)
        xlabel('Time (s)')
        ylabel('Response')
        legend({'Noiseless time series', 'Estimated broadband signal'}, 'Location', 'NorthWest')

end