function [broadband, params] = extractBroadband(signal, params)

% Compute time varying broadband envelope of a time series
% broadband = extractBroadband(x, srate, method, bands)
%
% Inputs
%   signal:  data (time x n) (n is number of channels or epochs)
%
%   params:  set of parameters defining how broadband is computed; see
%            s_simulateBroadband.m for descriptions

srate  = params.simulation.srate;
bands  = params.analysis.bands;

if isa(bands, 'cell')
    % Entire range for broadband
    band_rg  = bands{1}; 
    
    % Bin width 
    band_w   = bands{2}; 

    % All bins
    lb       = band_rg(1):band_w:band_rg(2)-band_w;
    ub       = lb+band_w;
    bands   = [lb; ub]';
end
    
% band pass filter each sub-band
bp  = zeros([size(signal) size(bands,1)]);
for ii = 1:size(bands,1)
    bp(:,:, ii) = butterpass_eeglabdata(signal,bands(ii,:),srate);
end

% if only one time series, then eliminate singleton dimension
%if size(bp, 2) == 1, bp = squeeze(bp); end

% which dimension represents the multiple bands?
%banddim = length(size(bp));
banddim = 3;

%%%% FUNCTIONS %%%%

% define whitening function
whiten = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));

% define averaging function
switch params.analysis.averagebandshow
    case 'geomean'
        mn = @(x,dim) geomean(x,dim);
        mnstr = 'geomean';
    case 'mean'
        mn = @(x,dim) mean(x,dim);
        mnstr = 'mean';
end

% define broadband computation
switch params.analysis.measure
    case 'amplitude'
        bb = @(x) abs(hilbert(x));
        bbstr = 'ampl';
    case 'power'
        bb = @(x) abs(hilbert(x)).^2;
        bbstr = 'power';
    case 'logpower'
        bb = @(x) log10(abs(hilbert(x)).^2);
        bbstr = 'logpower';
%     case 'logpower normalized' % dora method
%         bb = @(x) log10(abs(hilbert(x)).^2) - mean(log10(abs(hilbert(x)).^2));
%         bbstr = 'logpower norm';
    otherwise
        error('Unrecognized measure %s', params.analysis.measure)
end

%%%% COMPUTATIONS %%%%

% apply whitening?
switch params.analysis.whitenbands
    case 'yes'
        bp = whiten(bp);
        whitenstr = 'whiten';
    case 'no'
        whitenstr = [];
end

switch params.analysis.averagebandswhen
    case 'before hilbert'
        broadband = bb(mn(bp,banddim));
        methodstr = [bbstr ' ' mnstr ' ' whitenstr ' '];
    case 'after hilbert'   
        broadband = mn(bb(bp),banddim);
        methodstr = [mnstr ' ' bbstr ' ' whitenstr ' '];
end

% epoch the broadband time series
% nt = length(params.simulation.t);
% n  = params.simulation.n;
% broadband = reshape(broadband, nt, n);

params.analysis.methodstr = methodstr;

switch params.plot.on
    case 'yes'
        t = params.simulation.t/params.simulation.srate; 
        
        % plot individual band-pass filtered time series + envelopes
        singletrial_bp = squeeze(bp(:,1,:));
        singletrial_envelope = bb(singletrial_bp);
        if size(bp,3) > 1
            bandNumbers = [1  2 round(size(bands,1)/2) size(bands,1)];
        else
            bandNumbers = 1;
        end
        plotNames = {'Lowest band', 'Second band', 'Middle band', 'Highest band'};
        for ii = 1:length(bandNumbers)
            fH = figure; set(fH, 'Color', 'w'); hold on
            bpToPlot = singletrial_bp(:,bandNumbers(ii))/max(abs(singletrial_bp(:,bandNumbers(ii))));
            plot(t,bpToPlot, 'b');
            envToPlot = singletrial_envelope(:,bandNumbers(ii))/max(abs(singletrial_envelope(:,bandNumbers(ii))));
            plot(t,envToPlot, 'k','LineWidth', params.plot.lnwdth)
            set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl);
            title('Lowest band');
            legend({'Band-pass filtered time series', 'Hilbert envelope'},'Location', 'NorthWest')
            xlabel('Time (s)')
            ylabel('Amplitude')
            %set(gca, 'YLim', [-max(abs(singletrial_bp(:))) max(abs(singletrial_bp(:)))]);
            set(gca, 'YLim', [-1 1]);
            title([plotNames{ii} ' (' num2str(bands(bandNumbers(ii),1)) '-' num2str(bands(bandNumbers(ii),2)) 'Hz)']);
        end
        
        % plot the mean envelope across bands: single-trial 
        fH = figure; set(fH, 'Color', 'w'); hold on
        plot(t,mn(singletrial_envelope,2), 'LineWidth', params.plot.lnwdth);
        % compare with mean across trials
        test = bb(bp);
        testmn = mn(test,3);
        plot(t,mean(testmn,2), 'LineWidth', params.plot.lnwdth);
        set(gca, 'FontSize', params.plot.fontsz, 'XLim', params.plot.xl);
        set(gca, 'YLim', [0 max(abs(testmn(:)))]);
        title('Envelope averages across bands');
        legend({'Single-trial', 'Mean across trials'},'Location', 'NorthWest')
        xlabel('Time (s)')
        ylabel('Amplitude')
end
return

% switch method
%     case {1}
%         % -- Method 1: 'abs(hilbert(mean(whiten(bp(x)))))';
%         broadband = abs(hilbert(mean(whiten(bp),banddim)));
%         calc = 'abs(hilbert(mean(whiten(bp))))';
%     case {2}
%         % -- Method 2: 'abs(hilbert(mean(whiten(bp(x))))).^2';
%         broadband = abs(hilbert(mean(whiten(bp), banddim))).^2;
%         calc = 'abs(hilbert(mean(whiten(bp)))).^2';
%         
%     case {3}
%         % -- Method 3: geomean(abs(hilbert(whiten(bp(x))))
%         broadband = geomean(abs(hilbert(whiten(bp))),  banddim);
%         calc = 'geomean(abs(hilbert(whiten(bp)))';
%         
%     case {4}
%         % -- Method 4: geomean(abs(hilbert(whiten(bp(x)))).^2)
%         broadband = geomean(abs(hilbert(whiten(bp))).^2, banddim);
%         calc = 'geomean(abs(hilbert(whiten(bp))).^2)';
%         
%     case {5}
%         % -- Method 5: geomean(abs(hilbert(bp(x))).^2)
%         broadband = geomean(abs(hilbert(bp)).^2, banddim);
%         calc = 'geomean(abs(hilbert(bp)).^2)';
%         
%     case {6}
%         % -- Method 6: mean(abs(hilbert(whiten(bp(x)))).^2)
%         broadband = mean(abs(hilbert(whiten(bp))).^2, banddim);
%         calc = 'mean(abs(hilbert(whiten(bp))).^2)';
% end