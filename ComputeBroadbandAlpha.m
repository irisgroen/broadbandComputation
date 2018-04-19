% Compute Broadband using 3 methods, and compute alpha

% Sample data set
load('ECoG_sample_data_01.mat');
x = data;
t = -1000:2000; % time points for epoching
epochs = uint64(bsxfun(@plus,ones(210,1)*t,  onsets(1:2:end)));

% BB 1: 20 hz bins, avoiding 60 and 180 (max 210)
bands = [70 90;
    90 110;
    110 130;
    130 150;
    150 170;
    190 210];

bb = extractBroadband(x, srate, 4, bands);
epochedBB1 = bb(epochs);


% BB 2: 40 hz bins, avoiding 60 and 180 (max 230)
bands = [70 120;
    120 170;
    190 230];

bb = extractBroadband(x, srate, 4, bands);
epochedBB2 = bb(epochs);

% BB 3: single bin (70-170, no zscoring)
bands = [70 170];
bb = extractBroadband(x, srate, 2, bands);
epochedBB3 = bb(epochs);

% Alpha
bp = butterpass_eeglabdata(x,[10 20]);
alpha = abs(hilbert(bp));
epochedA = alpha(epochs);


figure, plot(t, zscore(mean(epochedBB1)), ...
    t, zscore(mean(epochedBB2)), ...
    t, zscore(mean(epochedBB3)), ...
    t, zscore(mean(epochedA)), 'k--', ...
    'LineWidth', 3)
legend('BB 1', 'BB 2', 'BB 3', 'Alpha')




