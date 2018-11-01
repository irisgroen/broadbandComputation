%% bandwidth comparison 

whichElectrodes = {'MO01'};%%, 'MO02'};%,'MO03','MO04'};
whichTrials = {'TWOPULSE-1'};%,'TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5'};

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
    set(gca, 'XLim', [0 0.75]);
end

