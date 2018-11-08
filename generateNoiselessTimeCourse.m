function [spikeRate, params] = generateNoiselessTimeCourse(params)

% Set defaults
if ~isfield(params.simulation, 'resp') || isempty(params.simulation.resp)
    params.simulation.resp = 'step';
end

t = params.simulation.t/params.simulation.srate; 
spontaneousRate = (1/params.simulation.srate)*100;
spikeRate = zeros(size(t))+spontaneousRate;

switch params.simulation.resp
    case 'boxcar'
        stimulus_idx = t>.4 & t <= .8; 
        baseline_idx = t>0 & t <=.4; 
        spikeRate(stimulus_idx) = 1;
        xl = [0 1];
        
    case 'pulse'
        spikeRate(t>.5 & t < .51) = 1;
        xl = [.4 0.6];
        
    case 'step'
        stimulus_idx = t>0.5 & t <= 1; 
        baseline_idx = t>0 & t <= 0.5; 
        spikeRate(stimulus_idx) = spikeRate(stimulus_idx) + 1;
        xl = [0 1]; %[0.4 0.6];
        
    case 'smallsteps'
        spikeRate(t>0.1)   = spikeRate(t>0.1) + 0.05; 
        spikeRate(t>0.25)  = spikeRate(t>0.25)+ 0.1; 
        spikeRate(t>0.5)   = spikeRate(t>0.5) + 0.3; 
        spikeRate(t>0.75)  = spikeRate(t>0.75)+ 0.55; 
        xl = [0 1];
        
	case 'bigsteps'
        spikeRate(t>0.1)   = spikeRate(t>0.1) + 1;
        spikeRate(t>0.25)  = spikeRate(t>0.25)+ 2;
        spikeRate(t>.5)    = spikeRate(t>0.5) + 3; 
        spikeRate(t>.75)   = spikeRate(t>0.75)+ 9; 
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
        spikeRate = sin(2*pi*t*f)+2+spontaneousRate;
        spikeRate(t<0) = spontaneousRate;
        xl = [.2 0.8];
        
    case 'noise'
        spikeRate = smooth(rand(size(t)), 50)+spontaneousRate;
        xl = [0 .5];
        
    case 'pred dn'
        tmp = load('predDn.mat', 'predDn');
        spikeRate = tmp.predDn+spontaneousRate;
        xl = [0 1];
        
    case 'level'
        if isfield(params.simulation.opt, 'level')
            level = params.simulation.opt.level;
        else, level = 1;
        end
        spikeRate(t>-1) = spikeRate(t>-1) + level;
        %spikeRate = spikeRate + level;
        xl = [-2 2]; %[0.4 0.6];
end

switch params.plot.on
    case 'yes'
        % Plot the noiseless time series
        figure();
        plot(t, spikeRate, 'b', 'LineWidth', params.plot.lnwdth);
        title('Noiseless time series (input)')
        set(gca, 'XLim', xl, 'FontSize', params.plot.fontsz)
        xlabel('Time (s)')
        ylabel('Spike Rate')
end

params.plot.xl = xl;
if exist('stimulus_idx','var');params.plot.stimulus_idx = stimulus_idx;end
if exist('baseline_idx','var');params.plot.baseline_idx = baseline_idx;end
end