% Dataset specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% Load preprocessed data

bandwidths = [2 5 10 20 40];

for jj = 1:length(bandwidths)
    fprintf('.');
    dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_epoched_bbmethod7_bandwidth%d', sub_label, ses_label, bandwidths(jj)));
    data(jj) = load(dataName);
end

%% bandwidth comparison 

whichElectrodes = {'MO02'};%%, 'MO02'};%,'MO03','MO04'};
whichTrials = {'TWOPULSE-1'};%,'TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5'};

smoothingLevelInMs = []; 
baselineType = 'all';

% 'out' contains the plotted time courses and SEs
for jj = 1:length(bandwidths)
    [out(jj)] = ecog_plotTrials(data(jj).trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs,baselineType);
    %set(gcf,'Name', [num2str(trials.bb_bands(1,2)-trials.bb_bands(1,1)) ' Hz bands']);
    a = get(gca);
    %a.Title.String = [whichElectrodes{:} '  :  ' data(bb_inx,jj).trials.bb_method '  :  ' num2str(data(bb_inx,jj).trials.bb_bands(1,2)-data(bb_inx,jj).trials.bb_bands(1,1)) ' Hz bands'];
    a.Title.String = [a.Title.String '  :  ' data(jj).trials.bb_method '  :  ' num2str(data(jj).trials.bb_bands(1,2)-data(jj).trials.bb_bands(1,1)) ' Hz bands'];
    %set(gca, 'XLim', [0 0.75]);
end

%% PLOT all conditions in one plot

colors = jet(length(bandwidths));
labels = [];
figure; hold on
for jj = 1:length(bandwidths)
    plot(out(jj).time,out(jj).broadband.(whichElectrodes{:}).mn, 'Color', colors(jj,:), 'LineWidth', 2)
    llim = out(jj).broadband.(whichElectrodes{:}).mn - out(jj).broadband.(whichElectrodes{:}).se;
    ulim = out(jj).broadband.(whichElectrodes{:}).mn + out(jj).broadband.(whichElectrodes{:}).se;
    h = ciplot(llim,ulim,out(jj).time,colors(jj,:), 0.25);
	h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    labels{jj} = ['bandwidth = ' num2str(bandwidths(jj))];
end
set(gca, 'XLim', [-0.25 0.75], 'FontSize', 18);
xlabel('Time(s)');
ylabel('Broadband response');
legend(labels);
% Add stim onset and zero lines
l = line([0 0], get(gca,'Ylim'),'LineStyle', ':', 'Color', 'k', 'LineWidth', 2);
l.Annotation.LegendInformation.IconDisplayStyle = 'off';
l = line([0.133333+0.016667 0.133333+0.016667],get(gca,'Ylim'),'LineStyle', ':', 'Color', 'k', 'LineWidth', 2);
l.Annotation.LegendInformation.IconDisplayStyle = 'off';
l = line([out(ii).time(1) out(ii).time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
l.Annotation.LegendInformation.IconDisplayStyle = 'off';
title(whichTrials);

