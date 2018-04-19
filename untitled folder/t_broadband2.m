%save ~/Desktop/twoChannPred.mat y
load ~/Desktop/twoChannPred.mat

srate0 = 1000;
nt = length(y);

% square pulse
y = zeros(size(y));
y(101:120) = 1;

% gaussian bump
y = zeros(size(y));
t = (1:length(y))./srate0;
t1 = 0.200; sigma = .030;
gau = exp(-(t-t1).^2/(2*sigma^2));
y = gau+.1;

n = 1000;
upsamplefactor = 5;
x = (1:nt)/srate0;
srate = srate0*upsamplefactor;
x2 = (1:nt*upsamplefactor)/(srate);
nt = length(x2);
rate = interp1(x,y, x2, 'spline');
rate = interp1(x,y, x2, 'linear', 'extrap');
%resp = bsxfun(@times,cumsum(randn(nt,n)),rate');
resp = bsxfun(@times,(randn(nt,n)),rate');
%spikes = poissrnd(ones(n,1)*rate, n, nt);

%resp = bsxfun(@times,randn(nt,n),rate');
data = resp(:);


%%

%% Load sample data
t = (1:size(data,1))/srate;
signal.raw = data; clear data;
f = (0:length(t)-1)/max(t);

%% Epochs
onsets = (0:n-2)*nt+1;
epochs = zeros(nt, length(onsets));
for ii = 1:length(onsets)
    epochs(:,ii) = onsets(ii)+(0:nt-1);
end

% Set first onset to t=0
t    = t-t(onsets(1));

% Define time for a single epoch
te  = t(onsets(1)+(0:nt-1))';

% Frequencies for a single epoch
fe = (0:length(te)-1)/max(te);


%% Plot time series
% figure(2); set(gcf, 'Color', 'w');
% 
% xl(1,:) = t([onsets(1) onsets(end)+nt]);
% xl(2,:) = t([onsets(1) onsets(4)]);
% for ii = 1:2
%     subplot(2,1,ii)
%     set(gca, 'FontSize', 18, 'XTick', t(onsets), 'XGrid', 'on', ...
%         'XTickLabel', round(t(onsets)))
%     hold on
%     plot(t, signal.raw);
%     xlabel('Time (s)'); ylabel('Amplitude (µV)')
%     xlim(xl(ii,:));
% end


% % Plot spectrum
% figure(3); set(gcf, 'Color', 'w');
% 
% set(gcf, 'Color', 'w'); set(gca, 'FontSize', 18), hold on
% plot(f, abs(fft(signal.raw))/length(t));
% plot(fe, abs(fft(mean(signal.raw(epochs),2)))/length(te), ...
%     'LineWidth', 4);
% set(gca, 'XScale', 'log', 'YScale', 'log')
% xlim([5 200])
% plot([60 60], get(gca, 'YLim'), 'k--')
% xlabel('Frequency (Hz)')
% ylabel('Amplitude (µV)')
% legend('Whole time series', 'Mean across epochs');


%% Define bands for subsequent band pass filtering
band_rg  = [100 200];
band_w   = 20;
lb       = band_rg(1):band_w:band_rg(2)-band_w;
ub       = lb+band_w;
bands   = [lb; ub]';
disp(bands)

%% Filter

% band pass filter: one wide band
signal.bp_single = butterpass_eeglabdata(signal.raw,band_rg,srate);

% band pass filter: several narrow bands
signal.bp_multi  = zeros(length(signal.raw),size(bands,1));
for ii = 1:size(bands,1)
    signal.bp_multi(:,ii) = butterpass_eeglabdata(signal.raw,bands(ii,:),srate);
    signal.bp_envelope(:,ii) = abs(hilbert(signal.bp_multi(:,ii)));
end

fH = figure(4); clf; set(fH, 'Color', 'w')
% make figure tall
sz = get(0, 'ScreenSize'); pos = get(fH, 'position');
set(fH, 'Position', [pos(1) sz(1) pos(3) sz(4)]);

% time for one epoch

for jj = 1:1
    clf
    idx = onsets(jj)+(0:nt-1);
    xl = [t(idx(1)) t(idx(end))];
    numplots = size(bands,1)+1;
    subplot(numplots,1,1)
    plot(t(idx), signal.bp_single(epochs), 'LineWidth', 1, 'Color', 'k');
    %plot(t(idx), signal.bp_single(idx), 'LineWidth', 1, 'Color', 'k');
    
    xlim(xl)
    title(sprintf('%d-%dHz', bands(1), bands(end)));
    colors = jet(numplots+1);
    for ii = 1:size(bands,1)
        subplot(numplots,1,ii+1)
        this_bp = signal.bp_multi(:,ii);
        plot(t(idx), this_bp(epochs), 'LineWidth', 1, 'Color', colors(ii,:));
        %plot(t(idx), this_bp(idx), 'LineWidth', 1, 'Color', colors(ii,:));
        hold on
        this_envelope = signal.bp_envelope(:,ii);
        %plot(t(idx), mean(this_envelope(epochs),2), 'LineWidth', 1, 'Color', 'k');
        plot(t(idx), mean(this_envelope(idx),2), 'LineWidth', 1, 'Color', 'k');
        
        %plot(t(idx), sum(signal.bp_envelope,2), 'LineWidth', 1, 'Color', 'g');
        
        xlim(xl)
        title(sprintf('%d-%d Hz', bands(ii,1), bands(ii,2)));
        
    end
    set(gcf, 'Name', sprintf('Epoch %d', jj));
    waitforbuttonpress
end


%% Compute time-varying broadband envelopes
methodnum = 0; str = [];
% 
% % -- Method 1: abs(hilbert) (no filtering) -----------------------------
% signal.bb{1} = abs(hilbert(signal.raw));
% str{1} = 'abs hilbert transform';
% 
% % --  Method 2: single bandpass, abs(hilbert) --------------------------
% signal.bb{2} = abs(hilbert(signal.bp_single));
% str{2} = 'single bandpass, abs(ht)';

% % -- Method 3: multiple bandpass, whiten, sum, abs(hilbert) ------------
% methodnum = methodnum+1;
% signal.bb{methodnum} = abs(hilbert(sum(1*(signal.bp_multi),2)));
% str{methodnum} = 'multiple bandpass, whiten, sum, abs(hilbert)';

% -- Method 4: multiple bandpass, abs(hilbert), whiten, sum ------------
methodnum = methodnum + 1;
signal.bb{methodnum} = sum(1*(abs(hilbert(signal.bp_multi))),2);
str{methodnum} = 'multiple bandpass, abs(ht), whiten, sum';

fH = figure(5);  clf; set(fH, 'Color', 'w')
set(gca, 'FontSize', 20)
hold on
numbb = length(signal.bb);



for ii = 1:numbb
    
    y = mean(signal.bb{ii}(epochs),2);
    plot(te, zscore(y), 'LineWidth', 4);
    
end

%
plot(t(idx), zscore(rate), 'k-', 'LineWidth', 3)

% blur the raw signal slightly
sigma = 0.010;
smoothing = 'gau'; % 'expo'

switch smoothing 
    case 'gau'
        kern   = exp(-(-.1:1/srate:.1).^2./(2*sigma^2));
        blurred = zscore(conv(rate, kern, 'same'));
    case 'expo'
        kern =  exp(-(0:1/srate:.1)./sigma);        
        blurred = zscore(conv(rate, kern, 'full'));
        some_shift = round(sigma*srate);
        blurred = blurred((1:length(rate))+some_shift);
end

plot(t(idx), blurred, 'g--', 'LineWidth', 3)

%plot(te, zscore(mean(signal.raw(epochs),2)), 'k', 'LineWidth', 2);
legend([str 'raw signal' 'blurred raw'], 'FontSize', 12, 'Location', 'Best')



return

