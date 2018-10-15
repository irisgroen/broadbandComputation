function [spikeRate, params] = generateNoiselessTimeCourse(params)

% Set defaults
if ~isfield(params.simulation, 'resp') || isempty(params.simulation.resp)
    params.simulation.resp = 'step';
end


t = params.simulation.t/params.simulation.srate; 
spontaneousRate = (1/params.simulation.srate)*10;
spikeRate = zeros(size(t))+spontaneousRate;


switch params.simulation.resp
    case 'boxcar'
        spikeRate(t>.2 & t < .7) = 1;
        xl = [0 1];
        
    case 'pulse'
        spikeRate(t>.5 & t < .51) = 1;
        xl = [.4 0.6];
        
    case 'step'
        spikeRate(t>0.5) = spikeRate(t>0.5) + 1;
        xl = [0 1];%[0.4 0.6];
        
    case 'steps'
        spikeRate(t>0.1)   = 1;
        spikeRate(t>0.25)  = 2;
        spikeRate(t>.5)    = 3;
        spikeRate(t>.75)   = 9;
        xl = [0 1];
        
    case 'bump'
        t1  = 0.500; sigma = .030;
        gau = exp(-(t-t1).^2/(2*sigma^2));
        spikeRate = gau+spontaneousRate; % + 1;
        xl = [.25 0.75];
        
    case 'square'
        if isfield(params.simulation.opt, 'f')
            f = params.simulation.opt.f;
        else, f = 3;
        end
        spikeRate = square(2*pi*t*f)+2+spontaneousRate;
        spikeRate(t<0) = spontaneousRate;
        xl = [.2 0.8];
        
    case 'sine'
        if isfield(params.simulation.opt, 'f')
            f = params.simulation.opt.f;
        else, f = 3;
        end
        spikeRate = sin(2*pi*t*f)+1+spontaneousRate;
        spikeRate(t<0) = spontaneousRate;
        xl = [.2 0.8];
        
    case 'noise'
        spikeRate = smooth(rand(size(t)), 50)+spontaneousRate;
        xl = [0 .5];
        
    case 'pred dn'
        tmp = load('predDn.mat', 'predDn');
        spikeRate = tmp.predDn+spontaneousRate;
        xl = [0 1];
        
end

switch params.plot.on
    case 'yes'
        % Plot the noiseless time series
        figure();
        plot(t, spikeRate, 'LineWidth', params.plot.lnwdth);
        title('Noiseless time series')
        set(gca, 'XLim', xl, 'FontSize', params.plot.fontsz)
        xlabel('Time (s)')
        ylabel('Spike Rate')
end

params.plot.xl = xl;

end