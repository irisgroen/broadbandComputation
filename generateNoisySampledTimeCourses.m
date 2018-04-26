
function [spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params)

% Generate multiple time series using poisson sampling at the time-varying rate determined above
t = params.t/params.srate; 

if ~isfield(params, 'seed') || isempty(params.seed)
    rng('shuffle')
else
    rng(params.seed,'twister'); 
end

spikeArrivals = poissrnd(repmat(spikeRate, [1 params.n]));

switch params.plot.on
    case 'yes'
        figure,
        plot(t, spikeArrivals(:,1), t, mean(spikeArrivals,2), 'LineWidth', params.plot.lnwdth);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)

        legend('Single noisy time series', 'Mean noisy time series')
        xlabel('Time (s)')
        ylabel('Spike Arrivals')
end

end
