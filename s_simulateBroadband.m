% s_simulateBroadband
%
% Simulate an LFP time series and then analyze it using our broadband
% extraction method. 
%
% SIMULATION: We simulate the LFP time series in 3 steps: 
% [1] There is an underlying, noiseless time series, which can be thought 
% of as the idealized rate of spike arrivals per neuron. 
% [2] We generate noisy samples scaled to this rate, which can be thought 
% of as a poisson-like nonlinearity defining the spike arrivals. 
% [3] We temporally integrate the noisy time series, which can be thought 
% of as the dendritic integration of the currents arising from spike arrivals.
%
% ANALYSIS: After generating multiple noisy time series with the
% identical underlying time-varying rate, we analyze the time series by
% extracting the broadband envelope in each trial using one of several
% algorithms. We then average the envelopes across trials computed by each
% algorithm, and compare these time-varying broadband envelopes to the
% noiseless time-varying rate used to seed the time series.

%% SIMULATION %%

% Choose the reponse profile and integration method; add noise (optional)
params = [];

% Set parameters for the noiseless, time-varying rate. 
params.resp        = 'steps';                 % response profile: choose from {'boxcar' 'steps' 'step' 'pulse' 'bump' 'square' 'sine' 'noise' 'pred dn'} ([default = step];
params.t           = (-1999.5:1999.5)';       % trial length: trials are -2 to 2 seconds, and later clipped to [0 1] to avoid edge artifacts
params.srate       = 1000;                    % sample rate (Hz) (Q shouldn't this go with the noisy sampling part? or would that be redundant)
% Set parameters for noisy samples
params.n           = 100;                     % number of repeated trials
params.seed        = 2;                       % use same number to compare simulations for same random generator of samples; leave empty to use new generator every time
% Set parameters for leaky integration
params.alpha       = 0.1;       % time constant for dendritic leakage
params.tau         = 0.0023;    % time constant for post-synaptic current
% Set parameters for noise
params.amplnoise   = 0.01;      % amplifier noise: scale factor of signal variance (if 0, no noise is added)
% Set parameters for plotting
params.plot.on     = 'yes';
params.plot.fontsz = 18; % font size
params.plot.lnwdth = 3;  % line width    

% [1] SIMULATE NOISELESS TIME SERIES
[spikeRate, params] = generateNoiselessTimeCourse(params);

% [2] GENERATE NOISY SAMPLES
[spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params);

% [3] TEMPORALLY INTEGRATE
[simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, params);

%% ANALYSIS %%

% This is an example analyses of one type of broadband computation

% [1] COMPUTE BROADBAND

% Define frequency bands and method for extracting broadband
params.bands       = {[70 170], 10}; % {[lower bound,  upper bound], window sz}
params.method      = 6;

[estimatedBroadband, params] = extractBroadband(simulatedSignal, params);

% [2] COMPARE WITH INPUT

[out] = evaluateBroadband(spikeRate, estimatedBroadband, params);

%% TO DO: COMPARISONS OF SIMULATIONS / ANALYSES %%

% Question: 
% How does the quality of broadband estimate vary using amplitude, power or log power estimates?
% Prediction: 
% Power best quality
% --> vary params.method

% Question: 
% How is temporal precision of broadband estimate affected by analysis parameters? 
% Prediction: 
% Time series containing sharp transients need small bands 
% --> vary params.resp / params.bands, where is the optimum?

% Question: 
% How is the quality of the broadband estimate affected by amplifier noise
% Prediction: 
% Quality decreases for high bands under conditions of high noise
% --> vary params.amplnoise / params.bands, where is the optimum?

% Implications for REAL data analysis:

% Given assumptions XY, e.g. amplitude vs power, will affect interpretation
% of selectivity (e.g. face response 2x or 4x as high as house response)
% How well can we capture actual transients in the data?
% Question: what are the BIG effects of choices made in analysis?
















