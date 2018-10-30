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

if isfield(params.simulation, 'stimulus_idx')
    stim_idx = params.simulation.stimulus_idx;
    base_idx = params.simulation.baseline_idx;
else 
    stim_idx = t > 0.1 & t < 0.4;
    base_idx = t > 0.6 & t < 0.9;
end

% Take the mean during level 0
y(1) = mean(mean(estimatedBroadband(stim_idx,:)));
x(1) = mean(spikeRate(stim_idx));

% Take the mean during level 1
y(2) = mean(mean(estimatedBroadband(base_idx,:)));
x(2) = mean(spikeRate(base_idx));

% Derive the slope & intercept
slope = diff(y)/diff(x);
intercept = y(2) - slope*(x(2));
calibrate = @(y) (y - intercept)./slope;

switch params.plot.on
    case 'yes'
        figure; title('calibration');
        scatter(x,y,200,'k', 'filled')
        hold on
        testy = linspace(y(1)-(y(2)),y(2)+y(2),20);
        testx = calibrate(testy);
        plot(testx,testy,'r-o', 'MarkerSize', 10, 'LineWidth', 2)
        xlabel('x (spikeRate)');
        ylabel('y (estimatedBroadband)')
        set(gca, 'FontSize', params.plot.fontsz)
        legend({'Calibrated points', 'Extrapolated points'}, 'Location', 'NorthWest');
        xl = get(gca, 'XLim');
        yl = get(gca, 'YLim');
        zerolines = plot(xl, [0 0], 'k--', [0 0], yl, 'k--');
        zerolines(1).Annotation.LegendInformation.IconDisplayStyle = 'off';     
        zerolines(2).Annotation.LegendInformation.IconDisplayStyle = 'off';     
end

% Add calibration function handle to params
params.analysis.calibrate = calibrate; 

end