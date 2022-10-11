% code for FC calculation based on coherence
clc; clear; close;
addpath(genpath('GrinstedWaveletCoherenceGus/'));

% load data signals and conduct coherence calculation
load('Projects\ROIsignals_all.mat')
connection = cell(size(signals, 1), 1);
parfor s = 1: size(signals, 1)
    % wtcMatrix
    %timeseries matrix:timpoint*ROI
    timeSeriesInput = squeeze(signals(s, :, :))';
    connection(s) = run_timeSeries2matFastCohi(timeSeriesInput, 2, [0.01 0.08], 190, 0, 'wavelet');     
    fprintf('%d is calculated!\n', s);
end
Coherence_connection = permute(cat(3, connection{:}), [3 1 2]);