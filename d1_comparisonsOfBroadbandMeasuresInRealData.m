% Dataset specs
projectName = 'visual';
sub_label   = 'som648'; 
ses_label   = 'nyuecog01';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';

% Output paths specs
dataDir = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label));

%% Load preprocessed data

bbmethods = [8 7 9];
bandwidths = [10 20 40];

for ii = 1:length(bbmethods)
    for jj = 1:length(bandwidths)
        fprintf('.');
        dataName = fullfile(dataDir, sprintf('sub-%s_ses-%s_epoched_bbmethod%d_bandwidth%d', sub_label, ses_label, bbmethods(ii), bandwidths(jj)));
        data(ii,jj) = load(dataName);
    end
    fprintf('\n');
end

%% broadband measure comparison

% Which stimulus conditions to plot? 
%whichElectrodes = {'IO01'};%, 'IO02','IO03','IO04'};
%whichTrials = {'SCENES','FACES', 'LETTERS'};%{'HRFPATTERN', 'SCENES','FACES', 'LETTERS'};

whichElectrodes = {'MO01'};%%, 'MO02'};%,'MO03','MO04'};
whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};

smoothingLevelInMs = [20]; 

% 'out' contains the plotted time courses and SEs across trials

% PLOT all conditions for each bb method separately
bw_inx = 1; %20hz bands
clear out;
for ii = 1:length(bbmethods)
    [out(ii)] = ecog_plotTrials(data(ii,bw_inx).trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs);
    %set(gcf,'Name', data(ii,bw_inx).trials.bb_method);
    a = get(gca);
    a.Title.String = [a.Title.String '  :  ' data(ii,bw_inx).trials.bb_method '  :  ' num2str(data(ii,bw_inx).trials.bb_bands(1,2)-data(ii,bw_inx).trials.bb_bands(1,1)) ' Hz bands'];
    set(gca, 'YLim', [-1 30]);
    set(gca, 'XLim', [0 0.75]);
end

% PLOT lowest and highest condition for all bb methods in one plot
figure; hold on
colors = copper(length(bbmethods));
colors(1,:) = [1 0 0];
for ii = 1:length(bbmethods)
    plot(out(ii).time,out(ii).broadband.(whichElectrodes{:}).mn(5,:),'-','Color', colors(ii,:),'LineWidth',2);
    plot(out(ii).time,out(ii).broadband.(whichElectrodes{:}).mn(1,:),':','Color', colors(ii,:),'LineWidth',2); 
end
ylim = get(gca,'YLim'); %ylim = [-1 30];
set(gca, 'YLim', ylim, 'XLim', [0 0.75], 'FontSize', 18);
% Add stim onset and zero lines
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([out(ii).time(1) out(ii).time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
legend({'amplitude CRF-5', 'amplitude CRF-1', 'power CRF-5', 'power CRF-1', 'logpower CRF-5', 'logpower CRF-1'});
xlabel('Time(s)')
ylabel(a.YLabel.String);
title('lowest and highest contrast for amplitude, power, logpower');

% PLOT lowest and highest condition for all bb methods in one plot
figure; hold on
colors = copper(length(bbmethods));
colors(1,:) = [1 0 0];
for ii = 1:length(bbmethods)
    plot(out(ii).time,out(ii).broadband.(whichElectrodes{:}).mn(1,:)/max(out(ii).broadband.(whichElectrodes{:}).mn(:)),':','Color', colors(ii,:),'LineWidth',2);
    plot(out(ii).time,out(ii).broadband.(whichElectrodes{:}).mn(5,:)/max(out(ii).broadband.(whichElectrodes{:}).mn(:)),'-','Color', colors(ii,:),'LineWidth',2); 
end
ylim = get(gca,'YLim'); %ylim = [-1 30];
set(gca, 'YLim', ylim, 'XLim', [0 0.75], 'FontSize', 18);
% Add stim onset and zero lines
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([out(ii).time(1) out(ii).time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
legend({'amplitude CRF-5', 'amplitude CRF-1', 'power CRF-5', 'power CRF-1', 'logpower CRF-5', 'logpower CRF-1'});
xlabel('Time(s)')
ylabel(a.YLabel.String);
title('lowest and highest contrast for amplitude, power, logpower');

% PLOT all conditions for power only with maxima highlighted
figure;hold on
bb_inx = 2; % power
colors = copper(length(whichTrials));
%colors(1,:) = [0.1 0.1 0.1];
for ii = 1:length(whichTrials)
    plot(out(bb_inx).time,out(bb_inx).broadband.(whichElectrodes{:}).mn(ii,:),'-','Color', colors(ii,:),'LineWidth',2);
    [y,x] = max(out(bb_inx).broadband.(whichElectrodes{:}).mn(ii,1:750)); 
    p = plot(out(bb_inx).time(x),y, 'Marker', '.', 'MarkerSize', 50, 'Color', colors(ii,:), 'LineStyle', 'none');
    p.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
ylim = [-1 30];
set(gca, 'YLim', ylim, 'XLim', [0 0.75], 'FontSize', 18);
% Add stim onset and zero lines
line([0 0], ylim,'LineStyle', ':', 'Color', 'k');
line([out(bb_inx).time(1) out(bb_inx).time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
legend(whichTrials);
xlabel('Time(s)')
ylabel(a.YLabel.String);
title('CRF 1-5 for power');

% PLOT maxima for all conditions, all bbmethods
maxima = []; maximaSE = [];
for ii = 1:length(whichTrials)
    for jj = 1:length(bbmethods)
        [maxima(ii,jj),x] = max(out(jj).broadband.(whichElectrodes{:}).mn(ii,1:750)); 
        maximaSE(ii,jj) = out(jj).broadband.(whichElectrodes{:}).se(ii,x);
    end
end

% compute the contrast of the stimuli
load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-som648_ses-nyuecog01_task-soc_run-01.mat')
CRFstimuli = stimulus.im_cell(find(contains(stimulus.categories, 'CRF')));
CRFrms = []; CRFp2p = [];
for ii = 1:length(CRFstimuli)
    CRFrms(ii) = rms(double(CRFstimuli{ii}(:))/255-0.5);
    CRFp2p(ii) = peak2peak(double(CRFstimuli{ii}(:))/255-0.5);
end

%% PLOT contrast response functions
figure;hold on
colors = copper(length(bbmethods));
colors(1,:) = [1 0 0];
for ii = 1:length(bbmethods)
    plot(CRFp2p,maxima(:,ii)/maxima(5,ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
    %plot(CRFp2p,maxima(:,ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
    %e = errorbar(CRFp2p,maxima(:,ii),maximaSE(:,ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
    e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlim([0 length(whichTrials)+1]);
end
%set(gca, 'YLim', [0 1.2], 'XLim', [0 1], 'FontSize', 18);
%set(gca, 'XTick', 1:length(whichTrials), 'XTickLabel', whichTrials);
title('contrast response functions');
xlabel('contrast level');
ylabel('maximum response');
legend('amplitude', 'power', 'logpower');
%set(gca, 'YScale','log')
set(gca, 'XScale','log')

