function [params] = calibrateResponseLevel(params)

% Apply the same parameters for simulation and analysis to a simple on-off
% time series, and use the resulting estimates of broadband to scale the
% calculated broadband response to the time series that is actually being
% tested.

% Copy parameters from main simulation
calibparams = params;

% Run a separate simulation on a simple step response profile
calibparams.simulation.resp        = 'step'; % has 0 for t > 0 & t < 0.5, has 1 for t > 0.5
calibparams.plot.on                = 'no'; 

[spikeRate, calibparams] = generateNoiselessTimeCourse(calibparams);
[spikeArrivals, calibparams] = generateNoisySampledTimeCourses(spikeRate, calibparams);
[simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, calibparams);
[estimatedBroadband, calibparams] = extractBroadband(simulatedSignal, calibparams);

t = calibparams.simulation.t/calibparams.simulation.srate; 

% Take the mean during level 0
idx = t > 0.1 & t < 0.4;
bbLevel0 = mean(mean(estimatedBroadband(idx),2));

% Take the mean during level 1
idx = t > 0.6 & t < 0.9;
bbLevel1 = mean(mean(estimatedBroadband(idx),2));

params.analysis.calibration = [bbLevel0 bbLevel1]; % two numbers, one for level 0, one for level 1

end