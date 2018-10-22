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
params.simulation.srate       = 1000;                    % sample rate (Hz) 

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
params.analysis.measure          = 'power';        % amplitude/power/logpower/

% PLOTTING parameters

params.plot.on = 'no';                                  % suppress plotting of simulations and each individual analysis
params.plot.fontsz = 18;                                 % font size
params.plot.lnwdth = 3;                                  % line width    

%% SIMULATE & ANALYZE

% COMPARE multiple temporal frequencies

params.simulation.resp        = 'sine';               
tempFrequencies               = [1:3:30];
windowSizes                   = {1, 5, 10, 25, 50};

bb = [];
stats = [];
for ii = 1:length(tempFrequencies)
    params.simulation.opt.f = tempFrequencies(ii);      
    [spikeRate, params] = generateNoiselessTimeCourse(params);
    [spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);
    [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);    
    for jj = 1:length(windowSizes)
        params.analysis.bands  = {[50 200], windowSizes{jj}};     % {[lower bound,  upper bound], window sz}   
        [bb{ii,jj}.out, bb{ii,jj}.params] = extractBroadband(simulatedSignal, params);
        [stats{ii,jj}, bb{ii,jj}.params] = evaluateBroadband(spikeRate, bb{ii,jj}.out, bb{ii,jj}.params);
    end
end

%% PLOT

rsqToPlot = [];
labels = [];
colors = jet(length(windowSizes));

for ii = 1:length(tempFrequencies)
    for jj = 1:length(windowSizes)
        rsqToPlot(ii,jj) = stats{ii,jj}.regress.rsq;
        labels{jj} = ['bandwidth = ' num2str(windowSizes{jj})];
    end
end
 
fH = figure;  set(fH, 'Color', 'w'); hold on;
for jj = 1:length(windowSizes)
    plot(tempFrequencies, rsqToPlot(:,jj), 'Color', colors(jj,:), 'Marker', 'o', 'LineWidth', params.plot.lnwdth);
end
set(gca, 'FontSize', params.plot.fontsz)

xlabel('Input frequency (Hz)')
ylabel('R2')
legend(labels, 'Location', 'NorthEast');
title('R2 with varying bandwidths for analysis');

t = params.simulation.t/params.simulation.srate;
% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;

avToPlot = [];
for ii = 1:length(tempFrequencies)
    for jj = 1:length(windowSizes)
        
        % Average across trials
        meanBroadband = mean(bb{ii,jj}.out,2);

        meanBroadbandCalibrated = bb{ii,jj}.params.analysis.calibrate(meanBroadband);

        meanBroadbandCalibrated = meanBroadbandCalibrated(idx);
        % Average across time
        avToPlot(ii,jj) = mean(meanBroadbandCalibrated);
        avToPlotSD(ii,jj) = std(meanBroadbandCalibrated,0,1);
    end
end
 
fH = figure;  set(fH, 'Color', 'w'); hold on;
for jj = 1:length(windowSizes)
    plot(tempFrequencies, avToPlot(:,jj), 'Color', colors(jj,:), 'Marker', 'o', 'LineWidth', params.plot.lnwdth);
end
set(gca, 'FontSize', params.plot.fontsz)

xlabel('Input frequency (Hz)')
ylabel('Mean broadband power')
legend(labels, 'Location', 'NorthEast');
title(['mean broadband ' bb{1,1}.params.analysis.measure ' with varying bandwidths']);
    