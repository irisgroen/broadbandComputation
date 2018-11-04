% q2_temporalPrecisionOfBroadbandAnalysis
%
% Question 2: temporal precision

% How is temporal precision of broadband estimate affected by analysis parameters? 
% Prediction: Time series containing sharp transients need wide bands (more time precision)
% Approach: vary params.resp and params.bands, where is the optimum?
clear all;
params = [];

% SIMULATION parameters

% Set parameters for the noiseless, time-varying rate 
params.simulation.resp        = 'sine';               % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.simulation.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.simulation.srate       = 1000;                    % sample rate (Hz) 

% Set parameters for noisy samples
params.simulation.nn          = 100;                     % number of neurons
params.simulation.ntrials     = 12;                      % number of repeated trials
params.simulation.seed        = 1;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;                     % time constant for dendritic leakage
params.simulation.tau         = 0.0023;                  % time constant for post-synaptic current

% Set parameters for noise
params.simulation.addnoise   = 0;%0.01;                   % additive noise (if 0, no noise is added)

% ANALYSIS parameters

params.analysis.averagebandshow  = 'geomean';             % geomean/mean
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
tempFrequencies               = 1:3:30;%[2:4:30];
windowSizes                   = {160};%{1, 5, 10, 20, 40, 80};

bb = [];
stats = [];
for ii = 1:length(tempFrequencies)
    params.simulation.opt.f = tempFrequencies(ii);      
    [spikeRate, params] = generateNoiselessTimeCourse(params);
    [spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);
    [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);    
    for jj = 1:length(windowSizes)
        params.analysis.bands  = {[40 200], windowSizes{jj}};     % {[lower bound,  upper bound], window sz}   
        [bb{ii,jj}.out, bb{ii,jj}.params] = extractBroadband(simulatedSignal, params);
        [stats{ii,jj}, bb{ii,jj}.params] = evaluateBroadband(spikeRate, bb{ii,jj}.out, bb{ii,jj}.params);
    end
end

%% PLOT RSQ

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

%% PLOT AMPLITUDE

t = params.simulation.t/params.simulation.srate;
% Clip time series to avoid edge artifacts
idx = t > 0 & t < 1;
f = 0:length(find(idx))-1; 

MnToPlot = []; 
VarToPlot = [];
AmpToPlot = [];

for ii = 1:length(tempFrequencies)
    for jj = 1:length(windowSizes)
        
        % Average across trials
        meanBroadband = mean(bb{ii,jj}.out,2);

        meanBroadbandCalibrated = bb{ii,jj}.params.analysis.calibrate(meanBroadband);

        meanBroadbandCalibrated = meanBroadbandCalibrated(idx);
        
        % Average across time
        MnToPlot(ii,jj) = mean(meanBroadbandCalibrated);
        VarToPlot(ii,jj) = var(meanBroadbandCalibrated,0,1);
        
        % Take spectrum, calculate amplitude
        spectrum = abs(fft(meanBroadbandCalibrated));
        AmpToPlot(ii,jj) = spectrum(f == tempFrequencies(ii));
    end
end
AmpToPlot = AmpToPlot/(length(f)/2); % scale to amplitude in signal

% amplitude
fH = figure;  set(fH, 'Color', 'w'); hold on;
p = plot(tempFrequencies, ones(length(tempFrequencies),1),'k--', 'LineWidth', 2);
for jj = 1:length(windowSizes)
    plot(tempFrequencies, AmpToPlot(:,jj), 'Color', colors(jj,:), 'Marker', '.', 'MarkerSize', 50, 'LineWidth', params.plot.lnwdth);
end
set(gca, 'FontSize', params.plot.fontsz)

xlabel('Input frequency (Hz)')
ylabel('Amplitude')
legend(['spikeRate amplitude' labels], 'Location', 'NorthEast');
title(['amplitude of broadband timeseries with varying bandwidths']);
  

% % amplitude
% f = 0:length(spikeRate(idx))-1; 
% figure, plot(f, abs(fft(spikeRate(idx))), 'LineWidth', params.plot.lnwdth); 
% figure, plot(f, abs(fft(meanBroadbandCalibrated)), 'LineWidth', params.plot.lnwdth); 

% %% average 
% fH = figure;  set(fH, 'Color', 'w'); hold on;
% p = plot(tempFrequencies, ones(length(tempFrequencies),1)*mean(spikeRate(idx)),'k--', 'LineWidth', 2);
% for jj = 1:length(windowSizes)
%     plot(tempFrequencies, MnToPlot(:,jj), 'Color', colors(jj,:), 'Marker', 'o', 'LineWidth', params.plot.lnwdth);
% end
% set(gca, 'FontSize', params.plot.fontsz)
% 
% xlabel('Input frequency (Hz)')
% ylabel('Mean broadband power')
% legend(['mean of spikeRate' labels], 'Location', 'NorthEast');
% title(['mean broadband ' bb{1,1}.params.analysis.measure ' with varying bandwidths']);
%     
% % variance
% fH = figure;  set(fH, 'Color', 'w'); hold on;
% p = plot(tempFrequencies, ones(length(tempFrequencies),1)*var(spikeRate(idx)),'k--', 'LineWidth', 2);
% for jj = 1:length(windowSizes)
%     plot(tempFrequencies, VarToPlot(:,jj), 'Color', colors(jj,:), 'Marker', 'o', 'LineWidth', params.plot.lnwdth);
% end
% set(gca, 'FontSize', params.plot.fontsz)
% 
% xlabel('Input frequency (Hz)')
% ylabel('Variance in broadband power')
% legend(['variance of spikeRate' labels], 'Location', 'NorthEast');
% title(['variance in broadband ' bb{1,1}.params.analysis.measure ' with varying bandwidths']);
%    

