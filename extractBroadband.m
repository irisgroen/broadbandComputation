function [broadband, params] = extractBroadband(signal, params)

% Compute time varying broadband envelope of a time series
% broadband = extractBroadband(x, srate, method, bands)
%
% Inputs
%   signal:      data (time x n) (n is number of channels or epochs)
%
%   params.srate:  sample rate (Hz) [default = 1000]
%
%   params.method: method for computing broadband, can be a number 1-6, or a
%            string. [default = 1]
%               string
%
%                 1 'abs(hilbert(mean(whiten(bp(x)))))';
%                 2 'abs(hilbert(mean(whiten(bp(x))))).^2';
%                 3 'geomean(abs(hilbert(whiten(bp(x)))))';
%                 4 'geomean(abs(hilbert(whiten(bp(x))).^2)';
%                 5 'geomean(abs(hilbert(bp)).^2'
%                 6 'mean(abs(hilbert(whiten(bp(x)))).^2)'
%   params.bands:  required for methods 2-6. Can be a matrix (number of bands x 2)
%                     or a cell array {[lb, ub], width}
%                       [default = {[60 200], 20]}


if ~isfield(params,'srate')  || isempty(params.srate),  srate = 1000; end
if ~isfield(params,'method') || isempty(params.method), method = 1;   end
if ~isfield(params,'bands')  || isempty(params.bands),  bands = {[60 200], 20}; end

srate  = params.srate;
method = params.method;
bands  = params.bands;

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

whiten = @(x) (x - mean(x(:)))./ diff(prctile(x, [.25 .75]));

switch method
    case {1}
        % -- Method 1: 'abs(hilbert(mean(whiten(bp(x)))))';
        broadband = abs(hilbert(mean(whiten(bp),banddim)));
        calc = 'abs(hilbert(mean(whiten(bp))))';
    case {2}
        % -- Method 2: 'abs(hilbert(mean(whiten(bp(x))))).^2';
        broadband = abs(hilbert(mean(whiten(bp), banddim))).^2;
        calc = 'abs(hilbert(mean(whiten(bp)))).^2';
        
    case {3}
        % -- Method 3: geomean(abs(hilbert(whiten(bp(x))))
        broadband = geomean(abs(hilbert(whiten(bp))),  banddim);
        calc = 'geomean(abs(hilbert(whiten(bp)))';
        
    case {4}
        % -- Method 4: geomean(abs(hilbert(whiten(bp(x)))).^2)
        broadband = geomean(abs(hilbert(whiten(bp))).^2, banddim);
        calc = 'geomean(abs(hilbert(whiten(bp))).^2)';
        
    case {5}
        % -- Method 5: geomean(abs(hilbert(bp(x))).^2)
        broadband = geomean(abs(hilbert(bp)).^2, banddim);
        calc = 'geomean(abs(hilbert(bp)).^2)';
        
    case {6}
        % -- Method 6: mean(abs(hilbert(whiten(bp(x)))).^2)
        broadband = mean(abs(hilbert(whiten(bp))).^2, banddim);
        calc = 'mean(abs(hilbert(whiten(bp))).^2)';
end

% epoch the broadband time series
nt = length(params.t);
n  = params.n;
broadband = reshape(broadband, nt, n);

params.methodstr = calc;

return