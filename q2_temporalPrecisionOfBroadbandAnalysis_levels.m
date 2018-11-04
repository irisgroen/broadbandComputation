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
params.simulation.nn          = 100;                     % number of neurons
params.simulation.ntrials     = 12;                      % number of repeated trials
params.simulation.seed        = 0;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;                     % time constant for dendritic leakage
params.simulation.tau         = 0.0023;                  % time constant for post-synaptic current

% Set parameters for noise
params.simulation.addnoise   = 0;%0.01;                   % additive noise (if 0, no noise is added)

% ANALYSIS parameters
params.analysis.averagebandshow  = 'geomean';          % geomean/mean
params.analysis.averagebandswhen = 'after hilbert';    % 'before hilbert'/'after hilbert'
params.analysis.whitenbands      = 'no';               % yes/no
params.analysis.measure          = 'power';
% PLOTTING parameters

params.plot.on = 'no';                                  % suppress plotting of simulations and each individual analysis
params.plot.fontsz = 18;                                % font size
params.plot.lnwdth = 3;                                 % line width    

%% compare different measures for different levels 

% DEFINE levels to be tested

inputLevels = 1;%((1:20)/20).^4; % range between 0 and 1
%inputLevels = ((1:20)/16).^4; % range between 0 and 2.44
%inputLevels = log((1:20/20)+2); 

% DEFINE bandwidths to be tested
windowSizes = {1, 5, 10, 20, 40, 80};

t = params.simulation.t/params.simulation.srate;

activityLevels   = nan(length(t), params.simulation.ntrials, length(inputLevels));
responseLevels   = nan(length(inputLevels), length(windowSizes));
responseLevelsSD = nan(length(inputLevels), length(windowSizes));

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

%% ANALYZE

% Clip time series to avoid edge artifacts
tidx = t > 0 & t < 1;
labels = [];
 meanBroadbandCalibrated = [];
for jj = 1:length(windowSizes)
    
    params.analysis.bands  = {[40 200], windowSizes{jj}};  
    [params] = calibrateResponseLevel(params);

    for ii = 1:length(inputLevels)
    
        simulatedSignal = activityLevels(:,:,ii);
    
        % Compute broadband
        [estimatedBroadband, params] = extractBroadband(simulatedSignal, params);
        
        % Average extracted broadband across trials
        %meanBroadband = mean(estimatedBroadband,2);

        % Scale broadband to calibrated response using calibrate function
        %meanBroadbandCalibrated(:,ii) = params.analysis.calibrate(meanBroadband);

        broadbandCalibrated = params.analysis.calibrate(estimatedBroadband);
        % Average level across time
        meanBroadbandCalibrated(:,ii) = mean(broadbandCalibrated(tidx,:),1);
        
        % Average and SD across trials
        responseLevels(ii,jj) = mean(meanBroadbandCalibrated(:,ii),1);
        responseLevelsSD(ii,jj) = std(meanBroadbandCalibrated(:,ii),0,1);

    end   
	labels{jj} = ['bandwidth = ' num2str(windowSizes{jj})];
    %figure;plot(t,squeeze(meanBroadbandCalibrated), 'LineWidth', 2);
    set(gca, 'FontSize', params.plot.fontsz)
    xlabel('Time (s)');
    ylabel('Broadband time course (calibrated)');
    title(labels{jj});
end

%% PLOT

fH = figure;  set(fH, 'Color', 'w'); hold on;
colors = jet(length(windowSizes));

for jj = 1:length(windowSizes)
    plot(inputLevels,responseLevels(:,jj), ...
        'LineWidth', 2, 'Marker', '.', 'MarkerSize', 25, 'Color', colors(jj,:))
    e = errorbar(inputLevels,responseLevels(:,jj), responseLevelsSD(:,jj), ... 
        'LineWidth', 2,'LineStyle', 'none', 'Color', colors(jj,:));
    e.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
legend(labels);
set(gca, 'XLim', [-0.1 max(inputLevels)+0.1*max(inputLevels)], 'FontSize', params.plot.fontsz)
xlabel('Input level (spike rate)');
ylabel('Broadband response level');
title('comparison of bandwidths');

fH = figure;  set(fH, 'Color', 'w'); hold on;
colors = jet(length(windowSizes));
level_inx = length(inputLevels);

p = plot(1:length(windowSizes), ones(length(windowSizes),1)* mean(spikeRate(tidx,level_inx),1),'k--', 'LineWidth', 2);
for jj = 1:length(windowSizes)
    plot(jj, responseLevels(level_inx,jj),'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 50, 'Color', colors(jj,:));
    e = errorbar(jj,responseLevels(level_inx,jj), responseLevelsSD(level_inx,jj), ... 
        'LineWidth', 2,'LineStyle', 'none', 'Color', colors(jj,:));
    e.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
set(gca, 'XTick', 1:length(windowSizes),'XTicklabel', windowSizes);
set(gca, 'XLim', [0 length(windowSizes)+1], 'FontSize', params.plot.fontsz);
%set(gca, 'YLim', [inputLevels(level_inx)-(0.5*inputLevels(level_inx)) inputLevels(level_inx)+(0.5*inputLevels(level_inx))]);
%legend(labels);
xlabel('BandWidths')
ylabel(['Estimated broadband for spike rate of ' num2str(inputLevels(level_inx))])
ylim([0.5 1.5])