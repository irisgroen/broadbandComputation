% q2_temporalPrecisionOfBroadbandAnalysis
%
% Question 2: temporal precision

% How is temporal precision of broadband estimate affected by analysis parameters? 
% Prediction: Time series containing sharp transients need wide bands (more time precision)
% Approach: vary params.resp and params.bands, where is the optimum?

params = [];

% SIMULATION parameters

% Set parameters for the noiseless, time-varying rate 
params.simulation.resp        = 'sine';               % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.simulation.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.simulation.srate       = 1000;                    % sample rate (Hz) (Q shouldn't this go with the noisy sampling part? or would that be redundant)
params.simulation.opt.f       = 10;                      % temporal frequency of response profile, applicable to sine wave or square wave

% Set parameters for noisy samples
params.simulation.n           = 100;                     % number of repeated trials
params.simulation.seed        = 1;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;                     % time constant for dendritic leakage
params.simulation.tau         = 0.0023;                  % time constant for post-synaptic current

% Set parameters for noise
params.simulation.amplnoise   = 0.01;                    % amplifier noise: scale factor of signal variance (if 0, no noise is added)

% ANALYSIS parameters

params.analysis.averagebandshow  = 'mean';             % geomean/mean
params.analysis.averagebandswhen = 'after hilbert';    % 'before hilbert'/'after hilbert'
params.analysis.whitenbands      = 'no';               % yes/no
params.analysis.measure          = 'power';    

% PLOTTING parameters
params.plot.on = 'no';                                  % suppress plotting of simulations and each individual analysis
params.plot.fontsz = 18;                                 % font size
params.plot.lnwdth = 3;                                  % line width    

%% SIMULATE

[spikeRate, params] = generateNoiselessTimeCourse(params);
[spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);
[simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);

%% ANALYZE: bandwidths

windowSizes = {1, 5, 10, 25, 50};
colors = jet(length(windowSizes));

bb = []; stats = [];
for ii = 1:length(windowSizes)
    params.analysis.bands  = {[50 200], windowSizes{ii}};     % {[lower bound,  upper bound], window sz}   
    [bb{ii}, params] = extractBroadband(simulatedSignal, params);
    [stats{ii}] = evaluateBroadband(spikeRate, bb{ii}, params); 
end

%% PLOT

fH = figure;  set(fH, 'Color', 'w'); hold on;
labels = [];

t = params.simulation.t/params.simulation.srate;
% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Plot spikeRate
spikeRateToPlot = spikeRate(idx) / norm(spikeRate(idx));
plot(t(idx), spikeRateToPlot, 'k:', 'LineWidth', params.plot.lnwdth)
labels{1} = 'idealized spike rate';

% Plot broadband
for ii = 1:length(windowSizes)
    
    meanBroadband = mean(bb{ii},2);
    
    % Subtract 'prestim' baseline
    baseline = meanBroadband(t > -1 & t < 0);
    meanBroadband = meanBroadband(idx) - mean(baseline);

    % Scale for plotting
    mnToPlot = meanBroadband / norm(meanBroadband);
    
    plot(t(idx), mnToPlot, 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth)
    set(gca, 'FontSize', params.plot.fontsz, 'XLim',  [0 1])
    xlabel('Time (s)')
    ylabel('Response')
    labels{ii+1} = ['bandwidth = ' num2str(windowSizes{ii}) ': r2 = ' num2str(round(stats{ii}.regress.rsq,2))];
end
legend(labels, 'Location', 'NorthWest');
title('temporal precision');

%% ANALYZE: lower bounds

lowerBounds = {20, 40, 60, 80, 100};
colors = jet(length(windowSizes));

bb = []; stats = [];
for ii = 1:length(lowerBounds)
    params.analysis.bands  = {[lowerBounds{ii} 200], 10};     % {[lower bound,  upper bound], window sz}   
    [bb{ii}, params] = extractBroadband(simulatedSignal, params);
    [stats{ii}] = evaluateBroadband(spikeRate, bb{ii}, params); 
end

%% PLOT

fH = figure;  set(fH, 'Color', 'w'); hold on;
labels = [];

t = params.simulation.t/params.simulation.srate;
% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Plot spikeRate
spikeRateToPlot = spikeRate(idx) / norm(spikeRate(idx));
plot(t(idx), spikeRateToPlot, 'k:', 'LineWidth', params.plot.lnwdth)
labels{1} = 'idealized spike rate';

% Plot broadband
for ii = 1:length(lowerBounds)
    
    meanBroadband = mean(bb{ii},2);
    
    % Subtract 'prestim' baseline
    baseline = meanBroadband(t > -1 & t < 0);
    meanBroadband = meanBroadband(idx) - mean(baseline);

    % Scale for plotting
    mnToPlot = meanBroadband / norm(meanBroadband);
    
    plot(t(idx), mnToPlot, 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth)
    set(gca, 'XLim', params.plot.xl);
    set(gca, 'FontSize', params.plot.fontsz, 'XLim',  [0 1])
    xlabel('Time (s)')
    ylabel('Response')
    labels{ii+1} = ['lowerbound = ' num2str(lowerBounds{ii}) ': r2 = ' num2str(round(stats{ii}.regress.rsq,2))];
end
legend(labels, 'Location', 'NorthWest');
title('temporal precision');

%% ANALYZE: lower bounds

upperBounds = {120, 140, 160, 180, 200};
colors = jet(length(windowSizes));

bb = []; stats = [];
for ii = 1:length(lowerBounds)
    params.analysis.bands  = {[50 upperBounds{ii}], 10};     % {[lower bound,  upper bound], window sz}   
    [bb{ii}, params] = extractBroadband(simulatedSignal, params);
    [stats{ii}] = evaluateBroadband(spikeRate, bb{ii}, params); 
end

%% PLOT
fH = figure;  set(fH, 'Color', 'w'); hold on;
labels = [];

t = params.simulation.t/params.simulation.srate;
% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Plot spikeRate
spikeRateToPlot = spikeRate(idx) / norm(spikeRate(idx));
plot(t(idx), spikeRateToPlot, 'k:', 'LineWidth', params.plot.lnwdth)
labels{1} = 'idealized spike rate';

% Plot broadband
for ii = 1:length(upperBounds)
    
    meanBroadband = mean(bb{ii},2);
    
    % Subtract 'prestim' baseline
    baseline = meanBroadband(t > -1 & t < 0);
    meanBroadband = meanBroadband(idx) - mean(baseline);

    % Scale for plotting
    mnToPlot = meanBroadband / norm(meanBroadband);
    
    plot(t(idx), mnToPlot, 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth)
    set(gca, 'XLim', params.plot.xl);
    set(gca, 'FontSize', params.plot.fontsz, 'XLim', [0 1])
    xlabel('Time (s)')
    ylabel('Response')
    labels{ii+1} = ['upperbound = ' num2str(upperBounds{ii}) ': r2 = ' num2str(round(stats{ii}.regress.rsq,2))];
end
legend(labels, 'Location', 'NorthWest');
title('temporal precision');