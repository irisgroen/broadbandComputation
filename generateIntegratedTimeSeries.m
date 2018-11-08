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
for ii = 1:params.simulation.ntrials
    Qconv = conv2(spikeArrivals(:,:,ii), psc, 'full')  ;
    Q(:,:,ii) = Qconv(1:size(spikeArrivals(:,:,ii), 1), :);
end

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

% Add noise?
if params.simulation.addnoise > 0
    noise = randn(size(I))*params.simulation.addnoise;
    I = I + noise;
end

% Average across neurons
signal = squeeze(mean(I,2));

% Plot?
switch params.plot.on
    case 'yes'
        
        signalMean = mean(signal,2);

        % plot single and sum of time series
        fH = figure;  set(fH, 'Color', 'w');
        plot(t, signal(:,1), t, signalMean, 'r', 'LineWidth', params.plot.lnwdth);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)

        legend('Single trial', 'Mean across trials', 'Location', 'NorthWest')
        xlabel('Time (s)')
        ylabel('Simulated Signal')
        title('LFP time series')

        if isfield(params.plot, 'stimulus_idx')
            % plot mean w different colors for baseline and stimulus
            fH = figure;  set(fH, 'Color', 'w'); hold on
            plot(t,signalMean, 'k', 'LineWidth', params.plot.lnwdth);
            plot(t(params.plot.stimulus_idx), signalMean(params.plot.stimulus_idx), 'r','LineWidth', params.plot.lnwdth); 
            set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl)
            legend('Baseline', 'Stimulus', 'Location', 'NorthWest')
            xlabel('Time (s)')
            ylabel('Simulated Signal')
            title('LFP time series')

            % plot frequency spectra for baseline and for stimulus
            fH = figure;  set(fH, 'Color', 'w');
            baseline = signal(params.plot.baseline_idx,:);
            stimulus = signal(params.plot.stimulus_idx,:);
            f = 0:length(baseline)-1; 
            plot(f, mean(abs(fft(baseline)),2),'k', f, mean(abs(fft(stimulus)),2), 'r', 'LineWidth', params.plot.lnwdth); 
            xlim([0 length(baseline)/2]); 
            %set(gca, 'XScale', 'log', 'YScale', 'log')
            set(gca, 'XLim', [0 200], 'YScale', 'log')
            set(gca, 'FontSize', params.plot.fontsz)
            legend('Baseline', 'Stimulus', 'Location', 'NorthWest')
            xlabel('Frequency(Hz)')
            ylabel('Power')
            title('Noiseless time series (input)')
            title('LFP spectra')

        end
end




end
