function [params] = calibrateResponseLevel(params)

% Apply the same parameters for simulation and analysis to a simple on-off
% time series, and use the resulting estimates of broadband to scale the
% calculated broadband response to the time series that is actually being
% tested.

% Copy parameters from main simulation
calibparams = params;

% Run a separate simulation on a simple step response profile
calibparams.simulation.resp        = 'step'; % has 0 for t > 0 & t < 0.5, has 1 for t > 0.5
calibparams.plot.on                = 'yes'; 

[spikeRate, calibparams] = generateNoiselessTimeCourse(calibparams);
[spikeArrivals, calibparams] = generateNoisySampledTimeCourses(spikeRate, calibparams);
[simulatedSignal] = generateIntegratedTimeSeries(spikeArrivals, calibparams);
[estimatedBroadband, calibparams] = extractBroadband(simulatedSignal, calibparams);

t = calibparams.simulation.t/calibparams.simulation.srate; 

% Take the mean during level 0
idx = t > 0.1 & t < 0.4;
y(1) = mean(mean(estimatedBroadband(idx),2));
x(1) = mean(spikeRate(idx));

% Take the mean during level 1
idx = t > 0.6 & t < 0.9;
y(2) = mean(mean(estimatedBroadband(idx),2));
x(2) = mean(spikeRate(idx));

% slope & intercept
slope = diff(y)/diff(x);
intercept = y(2) - slope*(x(2));
calibrate = @(y) (y - intercept)./slope;
params.analysis.calibration = calibrate; 

end