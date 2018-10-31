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

%% method comparison

%whichElectrodes = {'MO01'};%%, 'MO02'};%,'MO03','MO04'};
%whichTrials = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};

% Which stimulus conditions to plot? 
whichElectrodes = {'IO01'};%, 'IO02','IO03','IO04'};
whichTrials = {'SCENES','FACES', 'LETTERS'};%{'HRFPATTERN', 'SCENES','FACES', 'LETTERS'};

collapseTrialTypes = 'no'; 
smoothingLevelInMs = []; 

% 'out' contains the plotted time courses and SEs
bw_inx = 2; %20hz bands

for ii = 1:length(bbmethods)
    [out(ii)] = ecog_plotTrials(data(ii,bw_inx).trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs);
    %set(gcf,'Name', data(ii,bw_inx).trials.bb_method);
    a = get(gca);
    a.Title.String = [a.Title.String '  :  ' data(ii,bw_inx).trials.bb_method '  :  ' num2str(data(ii,bw_inx).trials.bb_bands(1,2)-data(ii,bw_inx).trials.bb_bands(1,1)) ' Hz bands'];
    set(gca, 'YLim', [-1 4]);
end

%% bandwidth comparison 

whichElectrodes = {'MO01'};%%, 'MO02'};%,'MO03','MO04'};
whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5'};

collapseTrialTypes = 'no'; 
smoothingLevelInMs = []; 

% 'out' contains the plotted time courses and SEs
bb_inx = 2; %power
for jj = 1:length(bandwidths)
    [out(jj)] = ecog_plotTrials(data(bb_inx,jj).trials, whichElectrodes, whichTrials, collapseTrialTypes, smoothingLevelInMs);
    %set(gcf,'Name', [num2str(trials.bb_bands(1,2)-trials.bb_bands(1,1)) ' Hz bands']);
    a = get(gca);
    %a.Title.String = [whichElectrodes{:} '  :  ' data(bb_inx,jj).trials.bb_method '  :  ' num2str(data(bb_inx,jj).trials.bb_bands(1,2)-data(bb_inx,jj).trials.bb_bands(1,1)) ' Hz bands'];
    a.Title.String = [a.Title.String '  :  ' data(bb_inx,jj).trials.bb_method '  :  ' num2str(data(bb_inx,jj).trials.bb_bands(1,2)-data(bb_inx,jj).trials.bb_bands(1,1)) ' Hz bands'];
end

