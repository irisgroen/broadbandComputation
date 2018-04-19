% Broadband tutorial
%
% Simulate a broadband time series by scaling white noise according to a
% time-varying rate. We do this for many repeated trials, and then compute
% the broadband envelope of the simulated time series on each trial. We
% then compare the average broadband envelope to the time-varying rate used
% to create the time series
%
% responses = {'step' 'pulse' 'bump' 'square' 'sine' 'noise'};
% for ii = 1:length(responses)
%   resp = responses{ii};
%   t_broadband;
% end


% Possible response profiles
responses = {'step' 'pulse' 'bump' 'square' 'sine' 'noise'};

% Default to 'step'. By putting the variable assignment within the if
% statement, we are able to run the script several times by command line
if ~exist('resp', 'var') || isempty(resp)
    resp = responses{1};
end

% noisetype
noisetype = 'poisson'; % white or poisson

% integration
integrationhow = 'leaky'; % 'cumsum' 'leaky' 'none'

srate  = 1000;  % sample rate
n      = 200;   % number of trials
t      = (-999.5:1999.5)'/srate;
nt     = length(t);  % number of time points

bands = {[80 200], 40}; % {[lb ub], bw}

% Generate noiseless time course
x      = zeros(size(t)); 

switch resp
    case 'pulse'        
        x(t>.5 & t < .52) = 1;
        xl = [.48 0.54];
        
    case 'step'
        x(t>0.5) = 1;
        x = x + 1;
        xl = [0.45 0.55];

    case 'bump'        
        t1  = 0.500; sigma = .030;
        gau = exp(-(t-t1).^2/(2*sigma^2));
        x   = gau + 1;
        xl = [.25 0.75];

   case 'square'
        x = square(2*pi*t*3)+2;
        xl = [.2 0.8];

   case 'sine'
        x = sin(2*pi*t*3)+2;
        xl = [.2 0.8];
        
    case 'noise'
        x = smooth(rand(size(t)), 50);
        xl = [0 .5];
end     

% Generate multiple noise samples scaled by x
switch noisetype
    case 'white'        
        noise = randn(nt, n);
        noise = bsxfun(@times, noise, x);
    case 'poisson'
        noise = poissrnd(repmat(x, [1 n]));
end

switch integrationhow
    case 'cumsum'
        y = cumsum(noise);
    case 'leaky'
        tau = 0.1;
        y = zeros(size(noise));
        y(1,:) = noise(1,:);
        for ii = 1:length(t)-1
            dy = (noise(ii,:) - y(ii,:)) / tau / srate;
            y(ii+1,:) = y(ii,:) + dy;
        end
    case 'none'
        y = noise;
end
%% Compute broadband using 4 different methods


str{1} = 'abs(hilbert(x))';
str{2} = 'abs(hilbert(bp(x)))';
str{3} = 'abs(hilbert(sum(whiten(bp(x)))))';
str{4} = 'sum(whiten(abs(hilbert(bp(x))))';

broadband{1} = extractBroadband(y(:), srate, 1);
broadband{2} = extractBroadband(y(:), srate, 2, bands{1});
broadband{3} = extractBroadband(y(:), srate, 3, bands);
broadband{4} = extractBroadband(y(:), srate, 4, bands);

% epoch the broadband time series
for ii = 1:4
    broadband{ii} = reshape(broadband{ii}, nt, n);    
end

%% Plot

fH = figure;  set(fH, 'Color', 'w');

% clip for plotting
idx = t > 0 & t < 1;


rescale = @(x) zscore(x);
for ii = 1:4       
    subplot(2,2,ii)
    set(gca, 'FontSize', 20)
    mn = mean(broadband{ii},2); 
    mn = mn(idx);
    
    % scale for plotting
    b = regress(x(idx), [mn ones(size(mn))]);
    mn = mn * b(1) + b(2);
    plot(t(idx), mn, t(idx), x(idx), 'k-', 'LineWidth', 4)
    title(str{ii})
    xlim(xl);
end
    



