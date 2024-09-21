function [zscored_signal] = trial_zscore(epoched_signal)
%trial_zscore Z-score an epoched signal sampled IID from a Gaussian across
%trials in dim 2
%   Detailed explanation goes here
num_trials = size(epoched_signal, 2);
mu = squeeze(mean(epoched_signal, 2));
centered_signal = epoched_signal - permute(mu, [1, 3, 2]);
sigma = std(epoched_signal, 0, 2);
zscored_signal = centered_signal ./ sigma;
end