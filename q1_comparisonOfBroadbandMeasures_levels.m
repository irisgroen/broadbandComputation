% q1_comparisonOfBroadbandMeasures
%
% Question 1: amplitude vs. power vs log power 
% How does the quality of broadband estimate vary using amplitude, power or log power estimates?
% Prediction: Power best quality
% Approach: vary params.method
clear all;

%% Set parameters 

params = [];

% SIMULATION parameters

% Set parameters for the noiseless, time-varying rate 
params.simulation.resp        = 'level';               % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.simulation.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.simulation.srate       = 1000;                    % sample rate (Hz) 
params.simulation.opt.f       = 10;                      % temporal frequency of response profile, applicable to sine wave or square wave
params.simulation.opt.level   = 1;                       % average level of response when using just one level

% Set parameters for noisy samples
params.simulation.n           = 100;                     % number of repeated trials
params.simulation.seed        = 1;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;                     % time constant for dendritic leakage
params.simulation.tau         = 0.0023;                  % time constant for post-synaptic current

% Set parameters for noise
params.simulation.amplnoise   = 0.01;%0.01;               % amplifier noise: scale factor of signal variance (if 0, no noise is added)

% ANALYSIS parameters

params.analysis.bands            = {[50 170], 20};     % {[lower bound,  upper bound], window sz}
params.analysis.averagebandshow  = 'mean';             % geomean/mean
params.analysis.averagebandswhen = 'after hilbert';    % 'before hilbert'/'after hilbert'
params.analysis.whitenbands      = 'no';               % yes/no

% PLOTTING parameters

params.plot.on = 'no';                                  % suppress plotting of simulations and each individual analysis
params.plot.fontsz = 18;                                % font size
params.plot.lnwdth = 3;                                 % line width    

%% compare different measures for different levels 

% DEFINE levels to be tested

inputLevels = ((1:20)/20).^4; % range between 0 and 1
inputLevels = ((1:20)/16).^4; % range between 0 and 2.44
%inputLevels = log((1:20/20)+2); 

powerMeasures = {'amplitude', 'power', 'logpower'};
t = params.simulation.t/params.simulation.srate;

activityLevels   = nan(length(t), params.simulation.n, length(inputLevels));
responseLevels   = nan(length(inputLevels), length(powerMeasures));
responseLevelsSD = nan(length(inputLevels), length(powerMeasures));

% SIMULATE
for ii = 1:length(inputLevels)
    
    params.simulation.opt.level = inputLevels(ii);
    
    [spikeRate(:,ii), params] = generateNoiselessTimeCourse(params);
    [spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate(:,ii), params);
    [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);
    
    activityLevels(:,:,ii) = simulatedSignal;
end

figure;plot(t,spikeRate, 'LineWidth', 2);
set(gca, 'FontSize', params.plot.fontsz)
xlabel('Time (s)');
ylabel('Spike rate');
title('Input levels');

figure;plot(t,squeeze(mean(activityLevels,2)), 'LineWidth', 2);
set(gca, 'FontSize', params.plot.fontsz)
xlabel('Time (s)');
ylabel('Integrated time series');
title('Simulated signal levels');

% ANALYZE

% Clip time series to avoid edge artifacts
tidx = t > 0 & t < 1;

for jj = 1:length(powerMeasures)
    
    params.analysis.measure = powerMeasures{jj};          
    % Get 'calibrated' responses for this powerMeasure
	[params] = calibrateResponseLevel(params);
  
    for ii = 1:length(inputLevels)
    
        simulatedSignal = activityLevels(:,:,ii);
    
        % Compute broadband
        [estimatedBroadband, params] = extractBroadband(simulatedSignal, params);
   
        % Average extracted broadband across trials
        meanBroadband = mean(estimatedBroadband,2);

        % Scale broadband to calibrated response using calibrate function
        meanBroadbandCalibrated(:,ii) = params.analysis.calibrate(meanBroadband);

        responseLevels(ii,jj) = mean(meanBroadbandCalibrated(tidx,ii),1);
        responseLevelsSD(ii,jj) = std(meanBroadbandCalibrated(tidx,ii),0,1);

    end   
    
    figure;plot(t,squeeze(meanBroadbandCalibrated), 'LineWidth', 2);
    set(gca, 'FontSize', params.plot.fontsz)
    xlabel('Time (s)');
    ylabel('Broadband time course (calibrated)');
    title(powerMeasures{jj});
end

%% PLOT

fH = figure;  set(fH, 'Color', 'w'); hold on;
colors = copper(length(powerMeasures));

for jj = 1:length(powerMeasures)
    plot(inputLevels,responseLevels(:,jj), ...
        'LineWidth', 2, 'Marker', '.', 'MarkerSize', 25, 'Color', colors(jj,:))
    e = errorbar(inputLevels,responseLevels(:,jj), responseLevelsSD(:,jj), ... 
        'LineWidth', 2,'LineStyle', 'none', 'Color', colors(jj,:));
    e.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
l = legend(powerMeasures, 'Location', 'NorthWest');
set(gca, 'XLim', [-0.1 max(inputLevels)+0.1*max(inputLevels)], 'FontSize', params.plot.fontsz)
xlabel('Input level (spike rate)');
ylabel('Broadband response level');
title('comparison of broadband measures');
