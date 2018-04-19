%% Broadband tutorial

%% Load sample data
load ECoG_sample_data_01
t = (1:size(data,1))/srate;
signal.raw = data; clear data;
f = (0:length(t)-1)/max(t);

%% Epoch the data
blank   = find(strcmp(stimnames, 'blank'));
targets = find(strcmp(stimnames, 'white noise') | ...
    strcmp(stimnames, 'pink noise') | ...
    strcmp(stimnames, 'brown noise'));
%targets = 1:7;
stim_onsets = round(onsets(ismember(stims, targets)));

all_onsets = round(onsets(stims ~= blank));
nt = median(diff(all_onsets));


epochs = zeros(nt, length(stim_onsets));
for ii = 1:length(stim_onsets)
    epochs(:,ii) = stim_onsets(ii)+(0:nt-1);
end

% Set first onset to t=0
t    = t-t(stim_onsets(1));

% Define time for a single epoch
t_e  = t(stim_onsets(1)+(0:nt-1))';

% Frequencies for a single epoch
f_e = (0:length(t_e)-1)/max(t_e);


%% Plot time series
figure(1); clf; pos = get(gcf, 'Position'); pos([3 4]) = [500 1000];
set(gcf, 'Color', 'w', 'Position', pos);


%  Plot the whole time series
subplot(4,1,1)
set(gca, 'FontSize', 12, 'XTick', t(stim_onsets), 'XGrid', 'on', ...
    'XTickLabel', ' '); hold on
plot(t, signal.raw);
xlabel('Time (s)'); ylabel('Amplitude (µV)')
xlim(t([stim_onsets(1) stim_onsets(end)+nt]));
title('Time series from whole experiment');

%  Plot the time series for a few epochs
subplot(4,1,2)
set(gca, 'FontSize', 12, 'XTick', t(stim_onsets), 'XGrid', 'on', ...
    'XTickLabel', round(t(stim_onsets))); hold on
plot(t, signal.raw, t, signal.raw);
xlabel('Time (s)'); ylabel('Amplitude (µV)')
xlim(t([stim_onsets(1) stim_onsets(4)]));
title('Time series from a few epochs');

% Plot spectrum 
subplot(4,1, 3:4)

set(gcf, 'Color', 'w'); set(gca, 'FontSize', 12), hold on
plot(f, abs(fft(signal.raw))/length(t));
plot(f_e, abs(fft(mean(signal.raw(epochs),2)))/length(t_e), ...
    'LineWidth', 4);
set(gca, 'XScale', 'log', 'YScale', 'log')
xlim([5 200])
plot([60 60], get(gca, 'YLim'), 'k--')
xlabel('Frequency (Hz)')
ylabel('Amplitude (µV)')
legend('Whole time series', 'Mean across epochs');
title('Amplitude spectrum');


%% Define bands for subsequent band pass filtering
band_rg  = [60 200];
band_w   = 20;
lb       = band_rg(1):band_w:band_rg(2)-band_w;
ub       = lb+band_w;
bands   = [lb; ub]';
disp(bands)

%% Filter

% band pass filter: one wide band
signal.bp_single = butterpass_eeglabdata(signal.raw,band_rg,srate);
signal.bp_envelope = abs(hilbert(signal.bp_single));

% band pass filter: several narrow bands
signal.bp_multi  = zeros(length(signal.raw),size(bands,1));
for ii = 1:size(bands,1)
    signal.bp_multi(:,ii) = butterpass_eeglabdata(signal.raw,bands(ii,:),srate);
end

% envelope of each band
signal.bp_envelopes = abs(hilbert(signal.bp_multi));

%% Plot the band pass filtered time series and envelope for one large band
fH = figure(3); clf; set(fH, 'Color', 'w')
% make figure wide
pos = get(fH, 'position'); pos([3 4]) = [1000 500];
set(fH, 'Position', pos);

% time for one epoch
idx = stim_onsets(1)+(0:nt-1);
xl = [t(idx(1)) t(idx(end))];
subplot(1,2,1)
plot(t(idx), signal.bp_single(idx), 'LineWidth', 1, 'Color', .6 * [1 1 1]); hold on
plot(t(idx), signal.bp_envelope(idx)*[1 -1], 'LineWidth', 1, 'Color', 'k');
yl = get(gca, 'YLim');
xlim(xl)
title(sprintf('%d-%dHz, one epoch', bands(1), bands(end)));

subplot(1,2,2)
plot(t(idx), signal.bp_single(epochs), 'LineWidth', 1, 'Color', .6 * [1 1 1]); hold on
plot(t(idx), mean(signal.bp_envelope(epochs),2)*[1 -1], 'LineWidth', 1, 'Color', 'k');

axis([xl yl])
title(sprintf('%d-%dHz, all epochs', bands(1), bands(end)));

%% Plot the band pass filtered time series and envelope for each band
fH = figure(4); clf; set(fH, 'Color', 'w')
pos = get(fH, 'position'); pos([3 4]) = [500 1000];
set(fH, 'Position', pos);

numplots = size(bands, 1);

colors = jet(numplots);
for ii = 1:numplots
    subplot(numplots,2,ii*2-1)
    this_ts = signal.bp_multi(:,ii);
    this_en = signal.bp_envelopes(:,ii);
    plot(t(idx), this_ts(idx), 'LineWidth', 1, 'Color', colors(ii,:)); hold on
    plot(t(idx), this_en(idx)*[1 -1], 'LineWidth', 1, 'Color', 'k');    
    xlim(xl); yl = get(gca, 'YLim');
    title(sprintf('%d-%d Hz', bands(ii,1), bands(ii,2)));
    
    subplot(numplots,2,ii*2)
    plot(t(idx), this_ts(epochs), 'LineWidth', 1, 'Color', colors(ii,:)); hold on
    plot(t(idx), mean(this_en(epochs),2)*[1 -1], 'LineWidth', 1, 'Color', 'k');
    
    axis([xl yl])
    title(sprintf('%d-%d Hz', bands(ii,1), bands(ii,2)));
    
end




%% Compute time-varying broadband envelopes

% -- Method 1: abs(hilbert) (no filtering) -----------------------------
signal.bb{1} = abs(hilbert(signal.raw));
str{1} = 'abs hilbert transform';

% --  Method 2: single bandpass, abs(hilbert) --------------------------
signal.bb{2} = abs(hilbert(signal.bp_single));
str{2} = 'single bandpass, abs(ht)';

% -- Method 3: multiple bandpass, whiten, sum, abs(hilbert) ------------
signal.bb{3} = abs(hilbert(sum(1*(signal.bp_multi),2)));
str{3} = 'multiple bandpass, whiten, sum, abs(hilbert)';

% -- Method 4: multiple bandpass, abs(hilbert), whiten, sum ------------
signal.bb{4} = sum(1*(abs(hilbert(signal.bp_multi))),2);
str{4} = 'multiple bandpass, abs(ht), whiten, sum';

fH = figure; set(fH, 'Color', 'w')
set(gca, 'FontSize', 20)
hold on
numbb = length(signal.bb);



for ii = 1:numbb
    
    y = mean(signal.bb{ii}(epochs),2);
    plot(t_e, zscore(y), 'LineWidth', 4);
    
end

plot(t_e, zscore(mean(signal.raw(epochs),2)), 'k', 'LineWidth', 2);
legend([str 'raw signal'], 'FontSize', 12, 'Location', 'Best')


return

