function [zscored_signal] = trial_zscore(epoched_signal)
%trial_zscore Z-score an epoched signal sampled IID from a Gaussian across
%trials in dim 3
%   Detailed explanation goes here
mu = squeeze(mean(epoched_signal, 3));
sigma = std(epoched_signal, 0, 3);
zscored_signal = (epoched_signal - mu) ./ sigma;
end