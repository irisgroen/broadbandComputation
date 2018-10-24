function [signal] = generateIntegratedTimeSeries(spikeArrivals, params)


% Simulate leaky integration in the dendrite to generate time-varying
% dendritic currents
alpha = params.simulation.alpha;
tau   = params.simulation.tau;

t = params.simulation.t/params.simulation.srate; 
dt = 1/params.simulation.srate;   % time step for simulations

% Post-synaptic current impulse response function. This current (multiplied
% by the synaptic weight) is initiated each time a spike arrives. It rises
% quickly and falls slowly. Normalize it to sum of 1.
psc = exp(-1/tau*(0:dt:.100))'; psc = psc / sum(psc);

% Generic post-synaptic current (Q)
Q = conv2(spikeArrivals, psc, 'full')  ;
Q = Q(1:size(spikeArrivals, 1), :);

% Weight the synaptic currents differently across synapes
synapticWeights = randn(1,size(Q,2));
synapticWeights = zscore(synapticWeights);
Qw = bsxfun(@times, Q, synapticWeights);

% Qwm = mean(Qw,2);
% Alternatively, don't model the post-synaptic currents

% Q = spikeArrivals;

% The leakage current (I) is proportional to the accumulated charge
I = zeros(size(Q));
%I(1,:) = spikeArrivals(1,:)+1;
for ii = 1:length(t)-1
    
    dIdt =  (-I(ii,:)/alpha +  Qw(ii,:) );
    
    dI   = dIdt * dt;
    
    I(ii+1,:) = I(ii,:) + dI;
end

I = I / alpha;

% Add amplifier output noise?
if params.simulation.amplnoise > 0
    %noise = randn(size(I))*std(I(:))*params.simulation.amplnoise;
    noise = randn(size(I))*params.simulation.amplnoise;
    signal = I + noise;
else
    signal = I;
end

% Plot?
switch params.plot.on
    case 'yes'
        
        figure,
        plot(t, signal(:,1), t, mean(signal,2), '--', 'LineWidth', params.plot.lnwdth);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)

        legend('Example of a single integrated time series', 'Mean across all integrated time series', 'Location', 'SouthEast')
        xlabel('Time (s)')
        ylabel('Simulated Signal')
end




end
