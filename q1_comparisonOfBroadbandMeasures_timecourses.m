% q1_comparisonOfBroadbandMeasures
%
% Question 1: amplitude vs. power vs log power 
% How does the quality of broadband estimate vary using amplitude, power or log power estimates?
% Prediction: Power best quality
% Approach: vary params.method

%% Set parameters 

params = [];

% SIMULATION parameters

% Set parameters for the noiseless, time-varying rate 
params.simulation.resp        = 'smallsteps';               % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.simulation.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.simulation.srate       = 1000;                    % sample rate (Hz) 
params.simulation.opt.f       = 10;                      % temporal frequency of response profile, applicable to sine wave or square wave

% Set parameters for noisy samples
params.simulation.n           = 100;                     % number of repeated trials
params.simulation.seed        = 1;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;                     % time constant for dendritic leakage
params.simulation.tau         = 0.0023;                  % time constant for post-synaptic current

% Set parameters for noise
params.simulation.amplnoise   = 0;% %0.01;                    % amplifier noise: scale factor of signal variance (if 0, no noise is added)

% ANALYSIS parameters

params.analysis.bands            = {[50 200], 10};     % {[lower bound,  upper bound], window sz}
params.analysis.averagebandshow  = 'mean';             % geomean/mean
params.analysis.averagebandswhen = 'after hilbert';    % 'before hilbert'/'after hilbert'
params.analysis.whitenbands      = 'no';               % yes/no

% PLOTTING parameters

params.plot.on = 'no';                                  % suppress plotting of simulations and each individual analysis
params.plot.fontsz = 18;                                % font size
params.plot.lnwdth = 3;                                 % line width    

%% SIMULATE

[spikeRate, params] = generateNoiselessTimeCourse(params);
[spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);
[simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);

%% ANALYZE

powerMeasures = {'amplitude', 'power', 'logpower'};
colors = copper(length(powerMeasures));

bb = []; stats = [];
for ii = 1:length(powerMeasures)
    params.analysis.measure = powerMeasures{ii};          
    [bb{ii}.out, params] = extractBroadband(simulatedSignal, params);
    [stats{ii}, bb{ii}.params] = evaluateBroadband(spikeRate, bb{ii}.out,params); 
end

%% PLOT

fH = figure;  set(fH, 'Color', 'w'); hold on;
labels = [];

t = bb{1}.params.simulation.t/bb{1}.params.simulation.srate;

% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

% Plot spikeRate
%spikeRateToPlot = spikeRate(idx) / norm(spikeRate(idx));
plot(t(idx), spikeRate(idx), 'k:', 'LineWidth', bb{1}.params.plot.lnwdth)
labels{1} = 'idealized spike rate';

for ii = 1:length(powerMeasures)
    
    % Get 'calibrated' responses to scale estimatedBroadband with
    [bb{ii}.params] = calibrateResponseLevel(bb{ii}.params);

    % Extracted broadband, averaged across trials
    meanBroadband = mean(bb{ii}.out,2);

    % Scale broadband to calibrated response
    meanBroadbandCalibrated = bb{ii}.params.analysis.calibrate(meanBroadband);

    plot(t(idx), meanBroadbandCalibrated(idx), 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth)
    labels{ii+1} = [powerMeasures{ii} ': r2 = ' num2str(round(stats{ii}.regress.rsq,2))];
end

set(gca, 'XLim', bb{1}.params.plot.xl);
set(gca, 'FontSize', bb{1}.params.plot.fontsz, 'XLim', [0 1])
xlabel('Time (s)')
ylabel('Response')
legend(labels, 'Location', 'NorthWest');
title('comparison of broadband measures');

%% to add: scatter plot of response levels vs. broadband level

