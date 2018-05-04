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
if size(bp, 2) == 1, bp = squeeze(bp); end

% which dimension represents the multiple bands?
banddim = length(size(bp));

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
    case 'logpower normalized' % dora method
        bb = @(x) log10(abs(hilbert(x)).^2) - mean(log10(abs(hilbert(x)).^2));
        bbstr = 'logpower norm';
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
    case 'before broadband'
        broadband = bb(mn(bp,banddim));
        methodstr = [bbstr ' ' mnstr ' ' whitenstr ' '];
    case 'after broadband'   
        broadband = mn(bb(bp),banddim);
        methodstr = [mnstr ' ' bbstr ' ' whitenstr ' '];
end

% epoch the broadband time series
nt = length(params.simulation.t);
n  = params.simulation.n;
broadband = reshape(broadband, nt, n);

params.analysis.methodstr = methodstr;

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