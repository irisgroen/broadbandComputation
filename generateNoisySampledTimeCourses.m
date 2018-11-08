
function [spikeArrivals, params] = generateNoisySampledTimeCourses(spikeRate, params)

% Generate multiple time series using poisson sampling at the time-varying rate determined above
t = params.simulation.t/params.simulation.srate; 

if ~isfield(params.simulation, 'seed') || isempty(params.simulation.seed)
    rng('shuffle')
else
    rng(params.simulation.seed,'twister'); 
end


spikeArrivals = poissrnd(repmat(spikeRate, [1 params.simulation.nn params.simulation.ntrials]));

switch params.plot.on
    case 'yes'
        figure, hold on,
        plot(t, spikeArrivals(:,1,1), 'c', 'LineWidth', params.plot.lnwdth);
        plot(t, mean(spikeArrivals(:,:,1),2), 'b', 'LineWidth', params.plot.lnwdth);
        plot(t, mean(mean(spikeArrivals,2),3), 'r', 'LineWidth', params.plot.lnwdth);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)

        legend('Single neuron, single trial', 'Single trial', 'Mean across trials', 'Location', 'NorthWest')
        xlabel('Time (s)')
        ylabel('Spike Arrivals')
        title('Spiking time series');
end

end
