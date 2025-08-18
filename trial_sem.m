function [normed_signal] = trial_sem(epoched_signal)
%trial_zscore Z-score an epoched signal sampled IID from a Gaussian across
%trials in dim 3
%   Detailed explanation goes here
sem = std(epoched_signal, 0, 3) ./ sqrt(size(epoched_signal, 3));
normed_signal = epoched_signal ./ sem;
end