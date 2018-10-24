% q3_effectsOfAmplifierNoise
%
% Question 3: Amplifier noise

% How is the quality of the broadband estimate affected by amplifier noise
% Prediction: Quality decreases for high bands under conditions of high noise
% Approach: vary params.amplnoise and params.bands, where is the optimum?

params = [];

% SIMULATION parameters

% Set parameters for the noiseless, time-varying rate 
params.simulation.resp        = 'smallsteps';               % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.simulation.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.simulation.srate       = 1000;                    % sample rate (Hz) 

% Set parameters for noisy samples
params.simulation.n           = 100;                     % number of repeated trials
params.simulation.seed        = 1;%[];                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time

% Set parameters for leaky integration
params.simulation.alpha       = 0.1;       % time constant for dendritic leakage
params.simulation.tau         = 0.0023;    % time constant for post-synaptic current

% ANALYSIS parameters

% Define frequency bands and method for extracting broadband
params.analysis.averagebandshow  = 'mean';             % geomean/mean
params.analysis.averagebandswhen = 'after hilbert';    % 'before hilbert'/'after hilbert'
params.analysis.whitenbands      = 'no';               % yes/no
params.analysis.measure          = 'power';            % amplitude/power/logpower

% PLOTTING parameters

params.plot.on     = 'no';
params.plot.fontsz = 12; % font size
params.plot.lnwdth = 3;  % line width    

%% SIMULATE 

[spikeRate, params] = generateNoiselessTimeCourse(params);
[spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);

%% SIMULATE % ANALYZE

% COMPARE RANGES
lowerBounds = [20 40 60 80 100]; % will add 100 for upper bound

bb = [];
stats = [];
for ii = 1:length(noiseLevels)
    params.simulation.amplnoise = noiseLevels(ii);      
    [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);
    for jj = 1:length(lowerBounds)
        params.analysis.bands = {[lowerBounds(jj) lowerBounds(jj)+100], 10};     
        %[bb{ii,jj}, params] = extractBroadband(simulatedSignal, params);
        %[stats{ii,jj}] = evaluateBroadband(spikeRate, bb{ii,jj}, params); 
        [bb{ii,jj}.out, bb{ii,jj}.params] = extractBroadband(simulatedSignal, params);
        [stats{ii,jj}, bb{ii,jj}.params] = evaluateBroadband(spikeRate, bb{ii,jj}.out, bb{ii,jj}.params); 
    end
end

%% PLOT

rsqToPlot = [];
labels = [];
for ii = 1:length(noiseLevels)
    for jj = 1:length(lowerBounds)
        rsqToPlot(ii,jj) = stats{ii,jj}.regress.rsq;
        labels{ii} = ['noise = ' num2str(noiseLevels(ii))];
    end
end
 
fH = figure;  set(fH, 'Color', 'w'); hold on;
for ii = 1:length(noiseLevels)
    plot(lowerBounds, rsqToPlot(ii,:), 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth);
end
set(gca, 'FontSize', params.plot.fontsz)
xnames = [];
for jj = 1:length(lowerBounds)
    xnames{jj} = [num2str(lowerBounds(jj)) '-' num2str(lowerBounds(jj)+100)];
end
set(gca, 'XTick', [min(lowerBounds):20:max(lowerBounds)], 'XTickLabel', xnames)
set(gca, 'XLim', [min(lowerBounds)-10 max(lowerBounds)+10], 'YLim', [0 1])
xlabel('Frequencies included in broadband analysis')
ylabel('R2')
legend(labels);    
title('Amplifier noise with varying frequency range for analysis');
set(gca, 'FontSize', params.plot.fontsz)

%% OLD 
% %% SIMULATE % ANALYZE
% 
% % GENERATE DIFFERENT NOISE REGIMES
% noiseLevels = [0 0.01 0.02 0.03 0.04 0.05];
% colors = parula(length(noiseLevels));
% 
% % COMPARE UPPER BOUNDS
% upperBounds = [100 120 140 160 180 200];
% 
% bb = [];
% stats = [];
% for ii = 1:length(noiseLevels)
%     params.simulation.amplnoise = noiseLevels(ii);      
%     [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);
%     for jj = 1:length(upperBounds)
%         params.analysis.bands = {[50 upperBounds(jj)], 10};     
%         [bb{ii,jj}, params] = extractBroadband(simulatedSignal, params);
%         [stats{ii,jj}] = evaluateBroadband(spikeRate, bb{ii,jj}, params); 
%     end
% end
% 
% %% PLOT
% 
% rsqToPlot = [];
% labels = [];
% for ii = 1:length(noiseLevels)
%     for jj = 1:length(upperBounds)
%         rsqToPlot(ii,jj) = stats{ii,jj}.regress.rsq;
%         labels{ii} = ['noise = ' num2str(noiseLevels(ii))];
%     end
% end
%  
% fH = figure;  set(fH, 'Color', 'w'); hold on;
% for ii = 1:length(noiseLevels)
%     plot(upperBounds, rsqToPlot(ii,:), 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth);
% end
% set(gca, 'FontSize', params.plot.fontsz)
% set(gca, 'XLim', [min(upperBounds)-10 max(upperBounds)+10], 'YLim', [0 1])
% xlabel('Upper bound')
% ylabel('R2')
% legend(labels);    
% title('Amplifier noise with varying upper bound');
% 
% %% SIMULATE % ANALYZE
% 
% % COMPARE LOWER BOUNDS
% lowerBounds = [20 40 60 80 100];
% 
% bb = [];
% stats = [];
% for ii = 1:length(noiseLevels)
%     params.simulation.amplnoise = noiseLevels(ii);      
%     [simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);
%     for jj = 1:length(lowerBounds)
%         params.analysis.bands = {[lowerBounds(jj) 150], 10};     
%         [bb{ii,jj}, params] = extractBroadband(simulatedSignal, params);
%         [stats{ii,jj}] = evaluateBroadband(spikeRate, bb{ii,jj}, params); 
%     end
% end
% 
% % PLOT
% 
% rsqToPlot = [];
% labels = [];
% for ii = 1:length(noiseLevels)
%     for jj = 1:length(lowerBounds)
%         rsqToPlot(ii,jj) = stats{ii,jj}.regress.rsq;
%         labels{ii} = ['noise = ' num2str(noiseLevels(ii))];
%     end
% end
%  
% fH = figure;  set(fH, 'Color', 'w'); hold on;
% for ii = 1:length(noiseLevels)
%     plot(lowerBounds, rsqToPlot(ii,:), 'Color', colors(ii,:), 'LineWidth', params.plot.lnwdth);
% end
% set(gca, 'FontSize', params.plot.fontsz)
% set(gca, 'XLim', [min(lowerBounds)-10 max(lowerBounds)+10], 'YLim', [0 1])
% xlabel('Lower bound')
% ylabel('R2')
% legend(labels);    
% title('Amplifier noise with varying lower bound');

